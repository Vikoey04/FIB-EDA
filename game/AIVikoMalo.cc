#include "Player.hh"

/**
 * Write the name of your player and save this file
 * with the same name and .cc extension.
 */
#define PLAYER_NAME VikoMalo

typedef vector<int> VI; //Vector d'enters


struct PLAYER_NAME : public Player {
    /**
     * Factory: returns a new instance of this class.
     * Do not modify this function.
     */
    static Player* factory () {
        return new PLAYER_NAME;
    }

    /**
    * Aquí es poden definir els tipus i atributs del teu jugador.
    */

    
    Dir BUSCA_PUNTO_CALIENTE(Pos pos) {
        map<Pos,Pos> visitades; //Set de posicions visitades
        queue<Pos> q;
        q.push(pos);
        visitades.insert(make_pair(pos,pos));
        Pos x;

        while (!q.empty()) {
            x = q.front();
            q.pop();
            
            for (int i = 0; i < 4; ++i) {
                //Si la pos està dins, si NO ES WALL i NO S'HA VISITAT
                if (pos_ok(x + Dir(i)) and cell(x + Dir(i)).type != WALL and visitades.find(x + Dir(i)) == visitades.end()) {
                    visitades.insert(make_pair(x + Dir(i), x));
                    q.push(x + Dir(i));

                    Cell c = cell(x + Dir(i));
                    if (c.type != GRASS) {
                        if (c.type == CITY and city_owner(c.city_id) != me()) {
                            map<Pos,Pos> :: iterator it = visitades.find(x + Dir(i));
                            x = it->second;
                            
                            while (x != pos) {
                                it = visitades.find(x);
                                x = it->second;
                            }
                            
                            if (it->first == (pos + BOTTOM)) return BOTTOM;
                            if (it->first == (pos + RIGHT)) return RIGHT;
                            if (it->first == (pos + TOP)) return TOP;
                            if (it->first == (pos + LEFT)) return LEFT;
                        }
                        else if (c.type == PATH and path_owner(c.path_id) != me()) {
                        
                            map<Pos,Pos> :: iterator it = visitades.find(x + Dir(i));
                            x = it->second;
                            
                            while (x != pos) {
                                it = visitades.find(x);
                                x = it->second;
                            }
                            
                            if (it->first == (pos + BOTTOM)) return BOTTOM;
                            if (it->first == (pos + RIGHT)) return RIGHT;
                            if (it->first == (pos + TOP)) return TOP;
                            if (it->first == (pos + LEFT)) return LEFT;
                        
                        } 
                    }
                    
                }
            }
            
            
        }
        return NONE;
        
    }
    

    /**
    * Mètode de joc, invocat un cop per ronda.
    */
    virtual void play () {
        VI U = my_units(me());
            int n = U.size();

            for (int i = 0; i < n; ++i) {
                int id = U[i];
                Unit u = unit(id); //Unitat
                //bool te_veins = false;

                for (int j = 0; j < 4; ++j) {
                    if (cell(u.pos+Dir(j)).unit_id != -1 and unit(cell(u.pos+Dir(j)).unit_id).player != me()) {
                        //Te un veí enemic en aquesta direccio, per tant l'ataquem:
                        move(id, Dir(j));
                    }
                    //Si no porta MASK i en té una al costat l'agafa
                    if (!u.mask) if (cell(u.pos+Dir(j)).mask) move(id, Dir(j));
                }
                
                //En aquest punt no té veïns, fem el BFS per trobar cami/ciutat més propera
                move(id, BUSCA_PUNTO_CALIENTE(u.pos));
            }
    }

};


/**
 * Do not modify the following line.
 */
RegisterPlayer(PLAYER_NAME);
