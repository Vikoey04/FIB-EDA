#include "Board.hh"
#include "Action.hh"


const char Board::uNDEF = gRASS;


void Board::generate_units () {
  pl_units_ = vector< vector<int> >(nb_players(), vector<int>(nb_units()));
  unit_ = vector<Unit>(nb_players() * nb_units());
  for (int id = 0, pl = 0; pl < nb_players(); ++pl) {
    for (int u = 0; u < nb_units(); ++u, ++id) {
      pl_units_[pl][u] = id;
      unit_[id].id = id;
      unit_[id].player = pl;
      unit_[id].health = initial_health();
      unit_[id].damage = (u < 3) ? u+1 : 0;
      unit_[id].turns = (u < 3) ? 1 : 0;
      unit_[id].immune = false;
      unit_[id].mask = false;
    }
    spawn(pl_units_[pl]);
  }
}


Board::Board (istream& is, int seed) {
  set_random_seed(seed);
  *static_cast<Settings*>(this) = Settings::read_settings(is);
  names_ = vector<string>(nb_players());
  read_generator_and_grid(is);

  round_ = 0;
  total_score_ = vector<int>   (nb_players(), 0);
  cpu_status_  = vector<double>(nb_players(), 0);

  city_owner_ = vector<int>(city_.size(), -1);
  path_owner_ = vector<int>(path_.size(), -1);
  generate_units();
  _my_assert(ok(), "Invariants are not satisfied.");
}


void Board::print_settings (ostream& os) const {
  // Should match the format of *.cnf files, except for the last line of board generation.
  os << version() << endl;
  os << "nb_players                  " << nb_players()                 << endl;
  os << "rows                        " << rows()                       << endl;
  os << "cols                        " << cols()                       << endl;
  os << "nb_rounds                   " << nb_rounds()                  << endl;
  os << "initial_health              " << initial_health()             << endl;
  os << "nb_units                    " << nb_units()                   << endl;
  os << "bonus_per_city_cell         " << bonus_per_city_cell()        << endl;
  os << "bonus_per_path_cell         " << bonus_per_path_cell()        << endl;
  os << "factor_connected_component  " << factor_connected_component() << endl;
  os << "infection_factor            " << infection_factor()           << endl;
  os << "mask_protection             " << mask_protection()            << endl;
}


void Board::print_names (ostream& os) const {
  os << "names         ";
  for (int pl = 0; pl < nb_players(); ++pl) os << ' ' << name(pl);
  os << endl;
}


void Board::print_state (ostream& os) const {

  // Should start with the same format of Info::read_grid.
  // Then other data describing the state.

  os << endl << endl;

  os << "   ";
  for (int j = 0; j < cols(); ++j)
    os << j / 10;
  os << endl;

  os << "   ";
  for (int j = 0; j < cols(); ++j)
    os << j % 10;
  os << endl;

  for (int i = 0; i < rows(); ++i) {
    os << i / 10 << i % 10 << " ";
    for (int j = 0; j < cols(); ++j) {
      const Cell& c = grid_[i][j];
      if (c.type == WALL) os << CellType2char(c.type);
      if (c.type == GRASS) {
	if (c.virus == 0) os << CellType2char(c.type);
	else os << char('a' + c.virus - 1);
      }
      if (c.type == PATH) {
	if (c.virus == 0) os << ',';
	else os << char('0' + c.virus - 1);
      }
      if (c.type == CITY) {
	if (c.virus == 0) os << ';';
	else os << char('A' + c.virus - 1);
      }
    }
    os << endl;
  }

  os << endl;
  os << "cities " << city_.size() << endl;
  for (int k = 0; k < int(city_.size()); ++k) {
    os << endl << city_[k].size() << endl;
    for (auto x: city_[k]) {
      os << x.i << " " << x.j << endl;
    }
  }

  os << endl;
  os << "paths " << path_.size() << endl;
  for (int k = 0; k < int(path_.size()); ++k) {
    os << endl
       << path_[k].first.first << " " << path_[k].first.second << " "
       << path_[k].second.size() << endl;
    for (auto x: path_[k].second) {
      os << x.i << " " << x.j << endl;
    }
  }

  os << endl;
  os << "masks " << masks_.size() << endl;
  for (int k = 0; k < int(masks_.size()); ++k) {
    os << masks_[k].i << " " << masks_[k].j << " " << endl;
  }

  os << endl;
  os << "round " << round() << endl;

  os << "total_score";
  for (auto ts : total_score_) os << " " << ts;
  os << endl;

  os << "status";
  for (auto st : cpu_status_) os << " " << st;
  os << endl;

  os << endl;
  os << "city_owners" << endl;
  for (int owner : city_owner_) os << " " << owner;
  os << endl;

  os << endl;
  os << "path_owners" << endl;
  for (int owner : path_owner_) os << " " << owner;
  os << endl;

  os << endl;
  os << "units" << endl;
  for (int id = 0; id < total_units(); ++id) {
    print_unit(unit(id), os);
    os << endl;
  }
  os << endl;
}


void Board::print_results () const {
  int max_score = 0;
  vector<int> v;
  for (int pl = 0; pl < nb_players(); ++pl) {
    cerr << "info: player " << name(pl)
         << " got score " << total_score(pl) << endl;
    if (total_score(pl) > max_score) {
      max_score = total_score(pl);
      v = vector<int>(1, pl);
    }
    else if (total_score(pl) == max_score) v.push_back(pl);
  }

  cerr << "info: player(s)";
  for (int pl : v) cerr << " " << name(pl);
  cerr << " got top score" << endl;
}


void Board::next(const vector<Action>& act, ostream& os) {

  _my_assert(ok(), "Invariants are not satisfied.");

  ++round_;

  int np = nb_players();
  int nu = total_units();

  // Chooses (at most) one command per unit.
  vector<bool> seen(nu, false);
  vector<Command> v;
  for (int pl = 0; pl < np; ++pl)
    for (const Command& m : act[pl].v_) {
      int id = m.id;
      Dir dir = m.dir;
      if (not unit_ok(id))
        cerr << "warning: id out of range : " << id << endl;
      else if (unit(id).player != pl)
        cerr << "warning: unit " << id << " of player " << unit(id).player
             << " not owned by " << pl << endl;
      else {
        // Here it is an assert because repetitions should have already been filtered out.
        _my_assert(not seen[id], "More than one command for the same unit.");
        seen[id] = true;
        if (not dir_ok(dir))
          cerr << "warning: direction not valid: " << dir << endl;
        else if (dir != NONE)
          v.push_back(Command(id, dir));
      }
    }

  // Executes commands using a random order.
  int num = v.size();
  vector<int> perm = random_permutation(num);
  vector<bool> killed(nu, false);
  vector<Command> commands_done;
  for (int i = 0; i < num; ++i) {
    Command m = v[perm[i]];
    if (not killed[m.id] and move(m.id, m.dir, killed))
      commands_done.push_back(m);
  }
  os << "commands" << endl;
  Action::print(commands_done, os);
  
  propagate(killed);

  // To ensure that executions of Game and SecGame are the same.
  for (int pl = 0; pl < np; ++pl)
    sort(pl_units_[pl].begin(), pl_units_[pl].end());

  vector<int> dead;
  for (int id = 0; id < nu; ++id)
    if (killed[id]) dead.push_back(id);

  spawn(dead);
  
  if (round_%5 == 0 and round_ < nb_rounds()) spawn_mask();

  compute_total_scores();

  _my_assert(ok(), "Invariants are not satisfied.");
}

// Spawns a mask in a GRASS cell without any units in it and without any mask on it
void Board::spawn_mask() {
	int at = 0;
	bool found = false;
	int i,j;
	while (at < MAX_ATTEMPTS and not found) {
		i = random(2, rows() - 3);
		j = random(2, cols() - 3);
		Cell& c = grid_[i][j];
		if (c.type == GRASS and c.unit_id == -1 and not c.mask) found = true;
	}
	if (found) {
		grid_[i][j].mask = true;
		masks_.push_back(Pos(i, j));
		return;
	}
	
	// If else fails, try exhaustively (this shouldn't happen)
	for (int i = 0; i < rows(); ++i) {
		for (int j = 0; j < cols(); ++j) {
			Cell& c = grid_[i][j];
			if (c.type == GRASS and c.unit_id == -1 and not c.mask) {
				c.mask = true;
				masks_.push_back(Pos(i, j));
				return;
			}
		}
	}
}

// Checks whether (i1,j1) and (i2,j2) are both indoors or outdoors
bool Board::same(int i1, int j1, int i2, int j2) {
	if (not pos_ok(i1,j1) or not pos_ok(i2,j2)) return false;
	Cell c1 = grid_[i1][j1];
	Cell c2 = grid_[i2][j2];
	if (c1.type == WALL or c2.type == WALL) return false;
	if (c1.type == GRASS) return c2.type == GRASS;
	return c2.type != GRASS;
}

// Propagates the virus, and develops the infection in units
void Board::propagate(vector<bool>& killed) {
	
	// Every non-masked infected propagates the virus
	for (int id = 0; id < (int)unit_.size(); ++id) {
		if (not killed[id] and unit_[id].damage > 0 and not unit_[id].mask) {
			Pos p = unit_[id].pos;
			Cell& c = grid_[p.i][p.j];
			c.virus += 3;
			//if (c.type == CITY or c.type == PATH) c.virus = min(10, c.virus);
			//if (c.type == GRASS) c.virus = min(4, c.virus);
		}
	}
	
	// Virus "travels" to adjacent cells
	vector<vector<Cell>> new_grid = grid_;
	for (int i = 0; i < rows(); ++i) {
		for (int j = 0; j < cols(); ++j) {
			Cell c = cell(i, j);
			if (c.type == WALL) continue;
			int vir = max(0,grid_[i][j].virus - 1);
			if (pos_ok(i+1,j) and same(i,j,i+1,j)) vir = max(vir, grid_[i+1][j].virus - 1);
			if (pos_ok(i-1,j) and same(i,j,i-1,j)) vir = max(vir, grid_[i-1][j].virus - 1);
			if (pos_ok(i,j+1) and same(i,j,i,j+1)) vir = max(vir, grid_[i][j+1].virus - 1);
			if (pos_ok(i,j-1) and same(i,j,i,j-1)) vir = max(vir, grid_[i][j-1].virus - 1);
			new_grid[i][j].virus = vir;
			if (new_grid[i][j].type == GRASS) new_grid[i][j].virus = min(vir, 4);
			else new_grid[i][j].virus = min(vir, 10);
		}
	}
	grid_ = new_grid;
	
	for (int id = 0; id < (int)unit_.size(); ++id) {
		Unit& u = unit_[id];
		if (not killed[id]) {
			// Infect susceptible units
			if (u.damage == 0 and not u.immune) {
				double p = grid_[u.pos.i][u.pos.j].virus/infection_factor();
				if (u.mask) p /= mask_protection();
				if (bernoulli(p)) {
					u.damage = random(2, 5);
					u.turns = 1;
				}
			}
			// Decide if units are no longer infected
			else if (not u.immune) {
				++u.turns;
				double p = 0.001*(u.turns*u.turns/16. + 1);
				if (bernoulli(p)) {
					u.damage = 0;
					u.immune = true;
				}
			}
		}
	}	
	
	// Deal damage to infected units
	for (int id = 0; id < (int)unit_.size(); ++id) {
		unit_[id].health -= unit_[id].damage;
		if (unit_[id].health < 0) {
			kill(id, random(0, 3), killed);
		}
	}	
}

void Board::spawn(const vector<int>& gen) {

  // Generate set of candidate positions for generation.
  set<Pos> cands;
  
  for (int i = 1; i + 1 < rows(); ++i) {
	  cands.insert(Pos(i, 1));
	  cands.insert(Pos(i, cols() - 2));
  }
  for (int j = 1; j + 1 < cols(); ++j) {
	  cands.insert(Pos(1, j));
	  cands.insert(Pos(rows() - 2, j));
  }
  
  // Regenerate killed units using valid candidate positions.
  for (int id : gen) {
    Pos pos = Pos(-1, -1);
    while (pos == Pos(-1, -1) and not cands.empty()) {
      int k = random(0, cands.size()-1);
      auto it = cands.begin();
      advance(it, k);
      if (valid_to_spawn(*it)) pos = *it;
      cands.erase(it);
    }
    if (pos == Pos(-1, -1)) // This should very very rarely happen.
      for (int i = 0; i < rows() and pos == Pos(-1, -1); ++i)
        for (int j = 0; j < cols() and pos == Pos(-1, -1); ++j)
          if (valid_to_spawn(Pos(i, j)))
            pos = Pos(i, j);
    _my_assert(pos != Pos(-1, -1), "Cannot find a cell to regenerate units");
    place(id, pos);
  }
}


bool Board::valid_to_spawn(Pos pos) {
  CellType t = cell(pos).type;
  if (t == WALL or t == CITY or t == PATH) return false;
  for (int d = 0; d <= NONE; ++d) {
    Pos pos2 = pos + Dir(d);
    if (pos_ok(pos2) and cell(pos2).unit_id != -1)
      return false;
  }
  return true;
}


void Board::place(int id, Pos p) {
  _my_assert(unit_ok(id), "Invalid identifier.");
  _my_assert( pos_ok( p), "Invalid position.");
  unit_[id].pos = p;
  grid_[p.i][p.j].unit_id = id;
}


void Board::kill(int id, int pl, vector<bool>& killed) {
  _my_assert(   unit_ok(id), "Invalid identifier.");
  _my_assert( player_ok(pl), "Invalid player.");
  _my_assert(not killed[id], "Cannot already be dead.");
  killed[id] = true;

  Unit& u = unit_[id];
  grid_[u.pos.i][u.pos.j].unit_id = -1;

  if (pl != u.player) {
    auto& o = pl_units_[u.player];
    auto it = find(o.begin(), o.end(), id);
    _my_assert(it != o.end(), "Cannot find id to kill.");
    swap(*it, *o.rbegin());
    o.pop_back();
    pl_units_[pl].push_back(id);
    u.player = pl;
  }
  u.pos    = Pos(-1, -1);
  u.health = initial_health();
  u.immune = false;
  u.mask = false;
  if (random(0, 4)) {
		u.damage = 0;
		u.turns = 0;
	}
	else {
		u.damage = random(2, 4),
		u.turns = 1;
	}
}


// id is a valid unit id, moved by its player, and d is a valid dir != NONE.
bool Board::move(int id, Dir dir, vector<bool>& killed) {
  _my_assert(unit_ok(id ), "Invalid identifier.");
  _my_assert( dir_ok(dir), "Invalid direction");
  _my_assert( dir != NONE, "Direction cannot be NONE");
  Unit& u = unit_[id];
  _my_assert(u.health >= 0, "Health cannot be negative.");
  Pos p1 = u.pos;
  _my_assert(pos_ok(p1), "Initial position in movement is not ok.");

  Cell& c1 = grid_[p1.i][p1.j];
  _my_assert(c1.type != WALL, "Initial position cannot be wall.");

  Pos p2 = p1 + dir;
  if (not pos_ok(p2)) return false;

  Cell& c2 = grid_[p2.i][p2.j];
  if (c2.type == WALL) return false;

  int id2 = c2.unit_id;
  if (id2 != -1) {
    _my_assert(unit_ok(id2), "Invalid identifier.");
    Unit& u2 = unit_[id2];
    _my_assert(u2.health >= 0, "Health cannot be negative.");
    if (u2.player == u.player) return false;
    int damage = random(25, 40);
    u2.health -= damage;
    if (u2.health < 0) kill(id2,  u.player, killed);
    else return false;
  }

  c1.unit_id = -1;
  c2.unit_id = id;
  u.pos = p2;
  if (c2.mask == true and u.mask == false) {
		u.mask = true;
		c2.mask = false;
		auto it = find(masks_.begin(), masks_.end(), p2);
		swap(*it, *masks_.rbegin());
    masks_.pop_back();
	}
  return true;
}


void Board::compute_scores_city_or_path(int bonus, const vector<Pos>& v, int& owner) {

  vector<int> sc(nb_players(), 0);
  for (Pos p : v) {
    int uid = cell(p).unit_id;
    if (uid != -1) ++sc[unit(uid).player];
  }
  int max_sc = 0;
  int max_pl = -1; // *Only* player with maximum score (-1 if more than one).
  for (int pl = 0; pl < nb_players(); ++pl) {
    if (sc[pl] > max_sc) {
      max_sc = sc[pl];
      max_pl = pl;
    }
    else if (sc[pl] == max_sc) max_pl = -1;
  }
  if (max_pl != -1) {     // Change of owner.
    total_score_[max_pl] += bonus * v.size();
    owner = max_pl;
  }
  else if (owner != -1) { // The owner is the same.
    total_score_[owner] += bonus * v.size();
  }
}


int size_of_connected_component_of(int u,
                                   const vector<vector<int>>& g,
                                   vector<bool>& mkd) {
  mkd[u] = true;
  int s = 1;
  for (int v : g[u])
    if (not mkd[v])
      s += size_of_connected_component_of(v, g, mkd);
  return s;
}


void Board::compute_scores_graph(int pl) {
  map<int, int> id;
  int sz = 0;
  for (int k = 0; k < int(city_.size()); ++k)
    if (city_owner_[k] == pl) {
      id[k] = sz;
      ++sz;
    }
  // There could be repeated edges, but there will be few.
  vector<vector<int>> g(sz);
  for (int k = 0; k < int(path_.size()); ++k) {
    int a = path_[k].first.first;
    int b = path_[k].first.second;
    if (path_owner_[k] == pl and
        city_owner_[a] == pl and
        city_owner_[b] == pl) {
      int ida = id[a];
      int idb = id[b];
      g[ida].push_back(idb);
      g[idb].push_back(ida);
    }
  }

  vector<bool> mkd(sz, false);
  for (int u = 0; u < sz; ++u)
    if (not mkd[u]) {
      int s = size_of_connected_component_of(u, g, mkd);
      _my_assert(0 <= s and s <= 25, "Unexpected size of connected component.");
      total_score_[pl] += factor_connected_component() * int(1 << s);
    }
}


void Board::compute_total_scores () {

  for (int k = 0; k < int(city_.size()); ++k)
    compute_scores_city_or_path(bonus_per_city_cell(), city_[k], city_owner_[k]);

  for (int k = 0; k < int(path_.size()); ++k)
    compute_scores_city_or_path(bonus_per_path_cell(), path_[k].second, path_owner_[k]);

  for (int pl = 0; pl < nb_players(); ++pl)
    compute_scores_graph(pl);
}


// ***************************************************************************

// S, E, N, W
const vector<int> Board::DIRI4 = { 1,  0, -1,  0};
const vector<int> Board::DIRJ4 = { 0,  1,  0, -1};

// SW, S, SE, E, NE, N, NW, W
const vector<int> Board::DIRI8 = {  1,  1,  1,  0, -1, -1, -1,  0};
const vector<int> Board::DIRJ8 = { -1,  0,  1,  1,  1,  0, -1, -1};

void Board::generator1 (const vector<int>& param) {

  _my_assert(param.empty(), "GENERATOR1 requires no parameter.");

  margin = min(rows(), cols()) / 5;

  bool found = false;
  while (not found) {
    m = vector<vector<char>>(rows(), vector<char>(cols(), uNDEF));
    city_.clear();
    path_.clear();

    fill_borders_with_walls();
    place_cities();
    place_old_cities();
    place_paths();
    place_walls();
    
    recode_temp();

    found = valid();
  }

  grid_ = vector< vector<Cell> >(rows(), vector<Cell>(cols()));
  for (int i = 0; i < rows(); ++i)
    for (int j = 0; j < cols(); ++j)
      grid_[i][j].type = char2CellType(m[i][j]);

  for (int k = 0; k < int(city_.size()); ++k)
    for (auto x: city_[k]) {
      _my_assert(grid_[x.i][x.j].type == CITY, "Mismatch with cities.");
      grid_[x.i][x.j].city_id = k;
    }

  for (int k = 0; k < int(path_.size()); ++k)
    for (auto x: path_[k].second) {
      _my_assert(grid_[x.i][x.j].type == PATH, "Mismatch with paths.");
      grid_[x.i][x.j].path_id = k;
    }
}
