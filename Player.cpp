Player.cpp
Type
C++
Size
5 KB (4,694 bytes)
Storage used
0 bytes
Location
version1
Owner
Theodor Panagiotakopoulos
Modified
Sep 15, 2019 by Theodor Panagiotakopoulos
Created
Sep 15, 2019
Add a description
Viewers can download

#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "vec.h"
#include "hmm_model.h"
namespace ducks
{
vector<vector<hmm_model>> hmm_model_list(COUNT_SPECIES,vector<hmm_model>());
int states = 13;
int observations = 8;

int * observations_of_bird(Bird bird){
    int array_size = bird.getSeqLength();
    int *obsv = new int[array_size];
    for(int i=0;i<bird.getSeqLength();i++){
        obsv[i] = bird.getObservation(i);
    }
    return obsv;
}

Player::Player()
{
    
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    std::ofstream myfile;
    myfile.open ("info.txt", std::ofstream::out | std::ofstream::app);
    myfile << "Round: " << pState.getRound() << ' '; // Round
    myfile << "Number of Birds: " << pState.getNumBirds() << "\n";
    for (int i=0;i<pState.getNumBirds();i++){
        Bird tmpBird = pState.getBird(i);
        myfile << "\tPrevious Movement of Bird: " << i << ", " << tmpBird.getSeqLength() << "\n\t\t";
        for (int j=0;j<tmpBird.getSeqLength();j++){
            myfile << tmpBird.getObservation(j) << ", ";
        }
        myfile << "\n";
    }
    myfile.close();
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */
    //
    // This line choose not to shoot
    return cDontShoot;

    //This line would predict that bird 0 will move right and shoot at it
    //return Action(0, MOVE_RIGHT);
}

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to guess the species of each bird.
     * This skeleton makes no guesses, better safe than sorry!
     */
    int birdnumber = pState.getNumBirds();
    std::vector<ESpecies> lGuesses(pState.getNumBirds());
    int *guess = new int[birdnumber];
    for(int i=0;i<pState.getNumBirds();i++){
        if(pState.getRound() == 0){
            //lGuesses[i] = SPECIES_RAVEN; // most common bird; 
            lGuesses[i] = static_cast<ESpecies>(rand()%(COUNT_SPECIES-1)); // random
            //std::cout << static_cast<ESpecies>(rand()%(COUNT_SPECIES-1)) << std::endl;
            guess[i] = rand()%(COUNT_SPECIES-1);
        }else{
            Bird bird = pState.getBird(i);
            //lGuesses[i] = SPECIES_UNKNOWN;
            int closest_model = 3;
            double highest_index = -1;
            for(int j=0;j<hmm_model_list.size();j++){
                for(int k=0;k<hmm_model_list[j].size();k++){
                    double x = hmm_model_list[j][k].validation_index(bird.getSeqLength(),ducks::observations_of_bird(bird));
                    //std::cerr << x << '|';
                    if(x>highest_index){
                        highest_index =x;
                       
                        closest_model = j;
                        guess[i] = j;
                    }
                }

            }
            std::cerr << "For bird " << i << " pos"<< highest_index << '/' << std::endl;;
            //guess[i] = closest_model;
            lGuesses[i] = static_cast<ESpecies>(closest_model);
            
            //std::cerr << lGuesses[i] << std::endl;

        }
    }
    std::cerr << "Proposed ";
    for(int i=0;i<birdnumber;i++){
        std:cerr << guess[i] + " ";
    }
    //delete guess;
    std::cerr << std::endl;
    return lGuesses;
}

void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
{
    /*
     * If you hit the bird you are trying to shoot, you will be notified through this function.
     */
    std::cerr << "HIT BIRD!!!" << std::endl;
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */
    std::cerr << "Real Species" ;
    for(int i=0;i<pSpecies.size();i++){
        std::cerr << pSpecies[i] << ' ';
        hmm_model HMM(states,observations);
        //HMM.A.print();
        //std::cerr << i << std::endl;
        Bird bird = pState.getBird(i);
        //std::cerr << pSpecies[i];
        //bird.~Bird();
        HMM.train(bird.getSeqLength(),ducks::observations_of_bird(bird));
        //HMM.p.print();
        //std::cerr << "list size: " << hmm_model_list[pSpecies[i]].size() << " bird:" << pSpecies[i] << " ";
        // distribute pi
        for(int i=0;i<HMM.p.cols;i++){
            HMM.p(0,i) = HMM.p(0,i)/HMM.p.cols + 0.5/HMM.p.cols; 
        }
        
        hmm_model_list[pSpecies[i]].push_back(HMM);
        //std::cerr << " inserted" << std::endl;
    }
    std::cerr << std::endl;
}



} /*namespace ducks*/
