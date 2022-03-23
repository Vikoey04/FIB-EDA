#include "Settings.hh"

Settings Settings::read_settings (istream& is) {
  // Should match the format of *.cnf files, except for the last line of board generation.
  Settings r;
  string s, v;

  // Version, compared part by part.
  istringstream vs(version());
  while (!vs.eof()) {
    is >> s;
    vs >> v;
    _my_assert(s == v, "Problems when reading.");
  };

  is >> s >> r.nb_players_;
  _my_assert(s == "nb_players", "Expected 'nb_players' while parsing.");
  _my_assert(r.nb_players_ == 4, "Wrong number of players.");

  is >> s >> r.rows_;
  _my_assert(s == "rows", "Expected 'rows' while parsing.");
  _my_assert(r.rows_ >= 20, "Wrong number of rows.");

  is >> s >> r.cols_;
  _my_assert(s == "cols", "Expected 'cols' while parsing.");
  _my_assert(r.cols_ >= 20, "Wrong number of columns.");

  is >> s >> r.nb_rounds_;
  _my_assert(s == "nb_rounds", "Expected 'nb_rounds' while parsing.");
  _my_assert(r.nb_rounds_ >= 1, "Wrong number of rounds.");

  is >> s >> r.initial_health_;
  _my_assert(s == "initial_health", "Expected 'initial_health' while parsing.");
  _my_assert(r.initial_health_ > 0, "Wrong initial health.");

  is >> s >> r.nb_units_;
  _my_assert(s == "nb_units", "Expected 'nb_units' while parsing.");
  _my_assert(r.nb_units_ >= 1, "Wrong number of units.");
  _my_assert(r.rows_ * r.cols_ >= 25 * r.nb_players_ * r.nb_units_, "Wrong parameters.");

  is >> s >> r.bonus_per_city_cell_;
  _my_assert(s == "bonus_per_city_cell", "Expected 'bonus_per_city_cell' while parsing.");
  _my_assert(r.bonus_per_city_cell_ >= 1, "Wrong bonus per city cell.");

  is >> s >> r.bonus_per_path_cell_;
  _my_assert(s == "bonus_per_path_cell", "Expected 'bonus_per_path_cell' while parsing.");
  _my_assert(r.bonus_per_path_cell_ >= 1, "Wrong bonus per path cell.");

  is >> s >> r.factor_connected_component_;
  _my_assert(s == "factor_connected_component", "Expected 'factor_connected_component' while parsing.");
  _my_assert(r.factor_connected_component_ >= 1, "Wrong factor for connected components.");
  
   is >> s >> r.infection_factor_;
  _my_assert(s == "infection_factor", "Expected 'infection_factor' while parsing.");
  _my_assert(r.factor_connected_component_ >= 1, "Wrong factor for infection.");
  
  is >> s >> r.mask_protection_;
  _my_assert(s == "mask_protection", "Expected 'mask_protection' while parsing.");
  _my_assert(r.factor_connected_component_ >= 1, "Wrong factor for mask protection.");
  
  _my_assert(r.rows_ == r.cols_, "Board should be square.");

  return r;
}
