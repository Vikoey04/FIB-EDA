#ifndef Board_hh
#define Board_hh


#include "Info.hh"
#include "Action.hh"
#include "Random.hh"


/*! \file
 * Contains the Board class, with all the game information,
 * plus the names of the players and a random generator.
 */


/**
 * Manages a random generator and the information of the board.
 */
class Board : public Info, public Random_generator {

  friend class Game;
  friend class SecGame;

  vector<string> names_;

  /**
   * Reads the generator method, and generates or reads the grid.
   */
  void read_generator_and_grid (istream& is) {
    string generator_;
    is >> generator_;
    if (generator_ == "FIXED") read_grid(is);
    else {
      vector<int> param;
      int x;
      while (is >> x) param.push_back(x);
      _my_assert(generator_ == "GENERATOR1", "Unknown grid generator.");
      generator1(param);
    }
  }


  /**
   * Prints some information of the unit.
   */
  inline static void print_unit (Unit u, ostream& os) {
    os << u.player << ' '
       << u.pos.i  << ' '
       << u.pos.j  << ' '
       << u.health << ' '
       << u.damage << ' '
       << u.turns << ' '
       << u.immune << ' '
       << u.mask << ' ';
  }

  bool valid_to_spawn(Pos pos);

  void place (int id, Pos p);
  
  void spawn(const vector<int>& gen);
  
  void kill(int id, int pl, vector<bool>& killed);
  
  void propagate(vector<bool>& killed);
  
  void spawn_mask();
  
  bool same(int i1, int j1, int i2, int j2);
  
  /**
   * Computes total scores of all players.
   */
  void compute_total_scores ();

  /**
   * Computes scores due to conquering a city/path with positions v.
   */
  void compute_scores_city_or_path(int bonus, const vector<Pos>& v, int& owner);

  /**
   * Computes scores due to the graph of player pl.
   */
  void compute_scores_graph(int pl);
  
  /**
   * Tries to apply a move. Returns true if it could.
   */
  bool move (int id, Dir dir, vector<bool>& killed);

  /**
   * Generate all the units of the board.
   */
  void generate_units ();


  /////////////////////// BEGIN BOARD GENERATION ///////////////////////

  static const char uNDEF;
  static const char tEMP = 't';

  /**
   * Generates a board.
   */
  void generator1 (const vector<int>& param);

  // Min and max number of *attempts* to place forests, etc.
  static const int MIN_NUM_CITIES      = 8;
  static const int MAX_NUM_CITIES      = 22;
  static const int MIN_NUM_WALLS       = 4;
  static const int MAX_NUM_WALLS       = 7;
  static const int MIN_NUM_OLD_CITIES  = 5;
  static const int MAX_NUM_OLD_CITIES  = 9;

  // Minimum distance ...
  static const int MIN_DISTANCE_OF_PATHS    =  2; // ... between paths.
  static const int MIN_DISTANCE_OF_CITIES   =  3; // ... between cities.
  static const int MIN_DISTANCE_OF_WALLS    =  2; // ... between walls.

  // Min and max horizontal and vertical sides of squares for cities.
  static const int MIN_CITY_HOR_SIDE  = 2;
  static const int MIN_CITY_VER_SIDE  = 2;
  static const int MAX_CITY_HOR_SIDE  = 6;
  static const int MAX_CITY_VER_SIDE  = 6;
  
  // Min and max sizes of sections of a wall
  static const int MIN_WALL_LENGTH = 5;
  static const int MAX_WALL_LENGTH = 10;

  // Max number of attempts before giving up placing an object.
  static const int MAX_ATTEMPTS = 100;
  
  // S, E, N, W
  static const vector<int> DIRI4;
  static const vector<int> DIRJ4;

  // SW, S, SE, E, NE, N, NW, W
  static const vector<int> DIRI8;
  static const vector<int> DIRJ8;
  
  int margin = 0;
  vector<vector<char>> m;
  vector<Pos> w;

  
  // Returns true with probability p.
  bool bernoulli(double p) {
    const int N = 10000;
    return random(1, N) <= int(N*p);
  }

  // Returns d such that central square is 
  // rows()/2 - d ... rows()/2 + d - 1 and
  // cols()/2 - d ... j <= cols()/2 + d - 1
  int half_side_of_central_square() {
    return ceil(sqrt(nb_units()));
  } 

  // Returns if (i, j) belongs to the central square.
  bool belongs_to_central_square(int i, int j) {
    int d = half_side_of_central_square();
    return
      i >= rows()/2 - d  and  i <= rows()/2 + d - 1 and
      j >= cols()/2 - d  and  j <= cols()/2 + d - 1;
  }

  // Returns whether there is a cell in the box of radius d
  // centered at (i, j) that contains character c.
  bool neighbour_one_of(const vector<vector<char>>& t, int i, int j, int d, char c) {
    for (int ii = max(0, i-d); ii <= min(rows()-1, i+d); ++ii)
      for (int jj = max(0, j-d); jj <= min(cols()-1, j+d); ++jj)
        if (t[ii][jj] == c)
          return true;
    return false;
  }

  // Returns whether there is a cell in the box of radius d
  // centered at (i, j) that contains an integer different from -1, a, b.
  bool neighbour_none_of(const vector<vector<int>>& t, int i, int j, int d, int a, int b) {
    for (int ii = max(0, i-d); ii <= min(rows()-1, i+d); ++ii)
      for (int jj = max(0, j-d); jj <= min(cols()-1, j+d); ++jj)
        if (t[ii][jj] != -1 and t[ii][jj] != a and t[ii][jj] != b)
          return true;
    return false;
  }

  // Returns the L^p distance between a and b.
  static double distance(const pair<double,double>& a, const pair<double,double>& b, double p = 2) {
    double first  = pow(abs(a.first  - b.first),  p);
    double second = pow(abs(a.second - b.second), p);
    return pow(first + second, 1/p);
  }

  // Returns if (i, j) is a position inside the map.
  bool inside(int i, int j) {
    return i >= 0 and i < rows() and j >= 0 and j < cols();
  }

  // Returns the number of elements in the matrix that are marked.
  static int size(const vector<vector<bool>>& mkd) {
    int sz = 0;
    for (int i = 0; i < int(mkd.size()); ++i)
      for (int j = 0; j < int(mkd[0].size()); ++j)
        sz += mkd[i][j];
    return sz;
  }

  // Returns a matrix where an area around position (i0, j0) is marked.
  template <class Prob>
  vector<vector<bool>> mark_area_around(int i0, int j0, Prob prob, bool allow_diags = true) {

    const vector<int>& diri = allow_diags ? DIRI8 : DIRI4;
    const vector<int>& dirj = allow_diags ? DIRJ8 : DIRJ4;

    vector<vector<bool>> mkd(rows(), vector<bool>(cols(), false));
    queue<Pos> q;
    q.push({i0, j0});
    mkd[i0][j0] = true;
    while (not q.empty()) {
      int i = q.front().i;
      int j = q.front().j;
      q.pop();
      for (int k = 0; k < int(diri.size()); ++k) {
        int ii = i + diri[k];
        int jj = j + dirj[k];
        if (inside(ii, jj) and not mkd[ii][jj] and prob(ii, jj)) {
          q.push({ii, jj});
          mkd[ii][jj] = true;
        }
      }
    }
    // Fill small gaps.
    bool changed = true;
    while (changed) {
      changed = false;
      for (int i = 2; i < rows()-2; ++i)
        for (int j = 2; j < cols()-2; ++j) {
          if (mkd[i-1][j] and not mkd[i][j] and mkd[i+1][j]) {
            mkd[i][j] = true;
            changed = true;
          }
          if (mkd[i][j-1] and not mkd[i][j] and mkd[i][j+1]) {
            mkd[i][j] = true;
            changed = true;
          }
          if (mkd[i-1][j] and not mkd[i][j] and not mkd[i+1][j] and mkd[i+2][j]) {
            mkd[i][j] = mkd[i+1][j] = true;
            changed = true;
          }
          if (mkd[i][j-1] and not mkd[i][j] and not mkd[i][j+1] and mkd[i][j+2]) {
            mkd[i][j] = mkd[i][j+1] = true;
            changed = true;
          }
        }
    }
    return mkd;
  }


  // Returns a vector of positions representing a curve starting at (i0, j0).
  template <class Prob>
  vector<Pos> curve_from(int i0, int j0, Prob prob, bool allow_diags = true) {

    const vector<int>& diri = allow_diags ? DIRI8 : DIRI4;
    const vector<int>& dirj = allow_diags ? DIRJ8 : DIRJ4;
    const int D = diri.size();

    vector<Pos> curve;
    vector<vector<bool>> mkd(rows(), vector<bool>(cols(), false));
    int i = i0;
    int j = j0;
    int k = random(0, D-1);
    while (true) {
      curve.push_back({i, j});
      mkd[i][j] = true;
      int s, ii, jj, kk;
      for (s = -2; s <= 1; ++s) {
        if (s < -1) kk = k + random(-1, 1); // First try random
        else        kk = k + s;             // then try exhaustively.
        // -1 <= kk - k <= 1 to get a smooth shape.
        kk = (D + kk) % D;
        _my_assert(0 <= kk and kk < D, "In curve generation.");
        ii = i + diri[kk];
        jj = j + dirj[kk];
        if (inside(ii, jj) and not mkd[ii][jj] and prob(i, j, ii, jj)) break;
      }
      if (s <= 1) { // Found a new point for the curve.
        i = ii;
        j = jj;
        k = kk;
      }
      else break; // Could not continue the curve, so stop.
    }
    return curve;
  }


  struct Prob1 {
    Board& b;
    int i0, j0;
    double prob;

    // Probability of getting true decreases as (i, j) moves away from (i0, j0).
    bool operator()(int i, int j) {
      return b.bernoulli(prob / distance({i, j}, {i0, j0}));
    }
  };


  struct Prob2 {
    Board& b;
    int i0, j0;
    double prob;

    // Probability of getting true decreases as (i, j) moves away from (i0, j0)
    // once (i, j) is far enough.
    bool operator()(int i, int j) {
      double p = b.random(20, 40)/10.; // p is a random real number between 2 and 4.
      double d = distance({i, j}, {i0, j0}, p);
      if (d <= min(b.rows(), b.cols())/2) d = 1;
      return b.bernoulli(prob / d);
    }
  };


  struct Prob3 {
    Board& b;
    int i0, j0;

    // True if (i, j) is closer to (i0, j0) than (ii, jj)
    bool operator()(int i, int j, int ii, int jj) {
      double  d = distance({i0, j0}, { i,  j});
      double dd = distance({i0, j0}, {ii, jj});
      return d < dd;
    }
  };


  struct Prob4 {
    Board& b;
    int i0, j0;

    // True if (ii, jj) is closer to (i0, j0) than (i, j)
    bool operator()(int i, int j, int ii, int jj) {
      double  d = distance({i0, j0}, { i,  j});
      double dd = distance({i0, j0}, {ii, jj});
      return d > dd;
    }
  };


  void fill_borders_with_walls() {
    for (int k = 0; k < rows(); ++k) m[k][0] = m[k][cols()-1] = wALL;
    for (int k = 0; k < cols(); ++k) m[0][k] = m[rows()-1][k] = wALL;
  }
  

  // Returns whether city with upper-left corner (i, j) and dimensions
  // (di, dj) is far enough of the sea, other cities or the central square.
  bool city_valid(int i, int j, int di, int dj) {
    const int D = MIN_DISTANCE_OF_CITIES;
    for (int ii = max(0, i-D); ii < min(rows(), i+di+D-1); ++ii)
      for (int jj = max(0, j-D); jj < min(cols(), j+dj+D-1); ++jj)
        if (m[ii][jj] == wALL or m[ii][jj] == cITY)
          return false;
    return true;
  }


  City place_city() {
    int i, j, di, dj;
    bool found = false;
    for (int k = 0; k < MAX_ATTEMPTS and not found; ++k) {
      i = random(0, rows()-1);
      j = random(0, cols()-1);
      di = random(MIN_CITY_VER_SIDE, MAX_CITY_VER_SIDE);
      dj = random(MIN_CITY_HOR_SIDE, MAX_CITY_HOR_SIDE);
      found = city_valid(i, j, di, dj);
    }
    City c;
    if (found) {
      for (int ii = i; ii < min(rows(), i+di); ++ii)
        for (int jj = j; jj < min(cols(), j+dj); ++jj)
          c.push_back({ii, jj});

      for (auto x : c)
        m[x.i][x.j] = cITY;
    }
    return c;
  }


  void place_cities() {
    int n_cities = random(MIN_NUM_CITIES, MAX_NUM_CITIES);
    for (int k = 0; k < n_cities; ++k) {
      auto c = place_city();
      if (not c.empty())
        city_.push_back(c);
    }
  }


  // Returns whether path p is too close to another existing path or
  // city (different from its leaving and arriving cities) or crosses
  // the central square.
  bool path_valid(const vector<Pos>& p,
               int a, int b,
               const vector<vector<int>>& mkd) {
    const int D = MIN_DISTANCE_OF_PATHS;
    for (auto x : p)
      if (neighbour_one_of(   m, x.i, x.j, D, pATH) or
					neighbour_one_of(   m, x.i, x.j, D, tEMP) or
          neighbour_none_of(mkd, x.i, x.j, D, a, b))
        return false;
    return true;
  }


  void place_paths() {

    vector<vector<int>> mkd(rows(), vector<int>(cols(), -1));
    int n_cities = city_.size();
    for (int k = 0; k < n_cities; ++k)
      for (const auto& x : city_[k])
        mkd[x.i][x.j] = k;

    int n_paths = 3*n_cities*n_cities;
    for (int k = 0; k < n_paths; ++k) {
      int a = random(0, city_.size()-1);
      int b = random(0, city_.size()-1);
      if (a != b) {                             // No auto-edges.
        int pa = random(0, city_[a].size()-1);
        int pb = random(0, city_[b].size()-1);
        int i1 = city_[a][pa].i;
        int j1 = city_[a][pa].j;
        int i2 = city_[b][pb].i;
        int j2 = city_[b][pb].j;
        auto c0 = curve_from(i1, j1, Prob4{*this, i2, j2}, false);
        if (c0.back() == Pos{i2, j2} and  // From (i1, j1) to (i2, j2).
            path_valid(c0, a, b, mkd)) {
          vector<Pos> c;
          for (auto x : c0)
            if (m[x.i][x.j] != cITY) { // Skin those cells from the
              m[x.i][x.j] = pATH;      // curve belonging to cities.
              c.push_back(x);
            }
          path_.push_back({{a, b}, c});
        }
      }
    }
  }
  
  // Returns whether old cities are far enough from everything
  bool old_city_valid(int i, int j, int di, int dj) {
    const int D = MIN_DISTANCE_OF_CITIES;
    for (int ii = max(0, i-D); ii < min(rows(), i+di+D-1); ++ii)
      for (int jj = max(0, j-D); jj < min(cols(), j+dj+D-1); ++jj)
        if (m[ii][jj] == wALL or m[ii][jj] == cITY or m[ii][jj] == pATH
         or m[ii][jj] == tEMP)
          return false;
    return true;
  }


  void place_old_city() {
    int i, j, di, dj;
    bool found = false;
    for (int k = 0; k < MAX_ATTEMPTS and not found; ++k) {
      i = random(0, rows()-1);
      j = random(0, cols()-1);
      di = random(MIN_CITY_VER_SIDE + 1, MAX_CITY_VER_SIDE);
      dj = random(MIN_CITY_HOR_SIDE + 1, MAX_CITY_HOR_SIDE);
      found = old_city_valid(i, j, di, dj);
    }
    if (found) {
      for (int ii = i; ii < min(rows(), i+di); ++ii) {
        for (int jj = j; jj < min(cols(), j+dj); ++jj) {
          m[ii][jj] = tEMP; // They will become walls later
          w.push_back(Pos(ii, jj));
				}
			}
    }
  }
  
  
  void place_old_cities() {
		int n_old_cities = random(MIN_NUM_OLD_CITIES, MAX_NUM_OLD_CITIES);
		for (int k = 0; k < n_old_cities; ++k) {
			place_old_city();
		}
	}
  
  
  bool wall_valid(vector<Pos>& c) {
		const int D = MIN_DISTANCE_OF_WALLS;
		for (auto x : c) {
			if (neighbour_one_of(m, x.i, x.j, D, wALL) or
					neighbour_one_of(m, x.i, x.j, D, pATH)) return false;
		}
		return true;
	}
	
	
	void place_walls() {
		int count = 0;
		if (w.size() == 0) return;
		for (int k = 0; k < MAX_ATTEMPTS and count < MAX_NUM_WALLS; ++k) {
			// Try to join two old cities or an old city with an existing one
			Pos p1, p2;
			p1 = w[random(0, w.size() - 1)];
			if (k%2 == 0) p2 = w[random(0,w.size() - 1)];
			else {
				int a = random(0, city_.size() - 1);
				p2 = city_[a][random(0, city_[a].size() - 1)];
			}
			int di = p1.i - p2.i;
			int dj = p1.j - p2.j;
			if (di*di + dj*dj > 2*MAX_CITY_HOR_SIDE*MAX_CITY_VER_SIDE) {
				auto c0 = curve_from(p2.i, p2.j, Prob4{*this, p1.i, p1.j}, false);
				double p = (double)random(55,75)/100.;
				vector<Pos> c;
				int l = 5;
				bool up = false;
				for (auto x : c0) {
					if (l == 0) {
						up = bernoulli(p); // Decide whether a segment of length l
						l = random(3,6);   // of the wall will be up or down
					}
					if (up) c.push_back(x);
					--l;
				}
				if (wall_valid(c)) {
					for (auto x : c) if(m[x.i][x.j] == gRASS) m[x.i][x.j] = wALL;
					++count;
				}
			}
		}
	}
	
	
	// Make old cities become walls
	void recode_temp() {
		for (int i = 0; i < rows(); ++i) 
			for (int j = 0; j < cols(); ++j) 
				if (m[i][j] == tEMP) m[i][j] = wALL;
	}


  // Mark all cells reachable by a unit from (i, j).
  void traversal(int i, int j, vector<vector<bool>>& mkd) {
    if (not mkd[i][j]) {
      mkd[i][j] = true;
      for (int k = 0; k < 4; ++k) {
        int ii = i + DIRI4[k];
        int jj = j + DIRJ4[k];
        if (inside(ii, jj) and m[ii][jj] != wALL)
          traversal(ii, jj, mkd);
      }
    }
  }


  // Returns whether a unit reach all cells without water.
  bool is_connected() {
    vector<vector<bool>> mkd(rows(), vector<bool>(cols(), false));
    traversal(rows()/2, cols()/2, mkd);

    for (int i = 0; i < rows(); ++i)
      for (int j = 0; j < cols(); ++j)
        if (m[i][j] != wALL and not mkd[i][j])
          return false;

    return true;
  }

  bool valid() {

    // Enough cities.
    if (city_.size() < MIN_NUM_CITIES) return false;

    // All cells without water should be reachable.
    if (not is_connected()) return false;

    return true;
  }
  
  /////////////////////// END BOARD GENERATION ///////////////////////

  
public:

  /**
   * Construct a board by reading information from a stream.
   */
  Board (istream& is, int seed);

  /**
   * Returns the name of a player.
   */
  inline string name (int player) const {
    _my_assert(player_ok(player), "Player is not ok.");
    return names_[player];
  }

  /**
   * Prints the board settings to a stream.
   */
  void print_settings (ostream& os) const;

  /**
   * Prints the name players to a stream.
   */
  void print_names (ostream& os) const;

  /**
   * Prints the state of the board to a stream.
   */
  void print_state (ostream& os) const;

  /**
   * Prints the results and the names of the winning players.
   */
  void print_results () const;

  /**
   * Computes the next board aplying the given actions to the current board.
   * It also prints to os the actual actions performed.
   */
  void next (const vector<Action>& act, ostream& os);

};

#endif
