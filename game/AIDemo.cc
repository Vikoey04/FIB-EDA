#include "Player.hh"


/**
 * Write the name of your player and save this file
 * with the same name and .cc extension.
 */
#define PLAYER_NAME Demo


// DISCLAIMER: The following Demo player is *not* meant to do anything
// sensible. It is provided just to illustrate how to use the API.
// Please use AINull.cc as a template for your player.

struct PLAYER_NAME : public Player {

  /**
   * Factory: returns a new instance of this class.
   * Do not modify this function.
   */
  static Player* factory () {
    return new PLAYER_NAME;
  }

  /**
   * Types and attributes for your player can be defined here.
   */

  typedef vector<int> VI;
  typedef vector<VI>  VVI;
  
  // Returns true if winning.
  bool winning() {
		for (int pl = 0; pl < nb_players(); ++pl)
      if (pl != me() and total_score(me()) <= total_score(pl))
        return false;
    return true;
  }
  
  // Moves all units of the player.
  void move_units() {
	  VI U = my_units(me()); // Get the id's of my units.
	  int n = U.size();
	  VI perm = random_permutation(n);
	  for (int i = 0; i < n; ++i) {
			
		  // We process the units in random order.
		  int id = U[perm[i]];
		  Unit u = unit(id);
		  
		  // With probability 1/4, we move at random.
		  if (random(0, 3)) move(id, Dir(random(0, DIR_SIZE - 1)));
		  
		  else if (u.damage > 0) { // If the unit is currently infected:
			  bool moved = false;
			  Pos p0 = u.pos;
			  for (int d = 0; d < 4 and not moved; ++d) {
				  Pos p1 = p0 + Dir(d);
				  if (pos_ok(p1)) {
						Cell c1 = cell(p1);
						if (c1.mask) { // If there's a mask nearby, go there.
							moved = true;
							move(id, Dir(d));
						}
					}
			  }
			  if (not moved) move(id, NONE);
		  }
		  else {
			  // Otherwise, do a bunch of (probably) stupid things.
			  // It's just to show that there are many possibilities.
			  if (u.immune) move(id, LEFT);
			  else if (cell(u.pos).type == CITY) {
					Cell c = cell(u.pos);
					if (city_owner(c.city_id) != me()) move(id, TOP);
					else if (c.virus >= 3) move(id, RIGHT);
					else move(id, BOTTOM);
				}
				else {
					if (cell(3, 4).type == WALL) move(id, NONE);
					else if (u.health < 30) move(id, BOTTOM);
					else if (cell(u.pos + TOP).unit_id != -1) {
						Unit u1 = unit(cell(u.pos + TOP).unit_id);
						if (u1.player == me()) move(id, RIGHT);
						else if (u1.turns >= 15) move(id, NONE);
						else if (u1.mask) move(id, LEFT);
						else move(id, TOP);
					}
					else if (random(0, 1)) move(id, Dir(random(0, 2)));
					// You can also use cerr to debug.
					// else cerr << u.pos.i << " " << u.pos.j << endl;
				}
		  }
	  }
  }
  
  /**
   * Play method, invoked once per each round.
   */
  virtual void play () {
    
    // If more than halfway through, do nothing.
    if (round() > nb_rounds()/2) return;
    
    // If winning, do nothing.
    if (winning()) return;

    // If nearly out of time, do nothing.
    double st = status(me());
    if (st >= 0.9) return;
    
    move_units();
  }

};


/**
 * Do not modify the following line.
 */
RegisterPlayer(PLAYER_NAME);

