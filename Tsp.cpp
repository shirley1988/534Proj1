#include <iostream>  // cout
#include <fstream>   // ifstream
#include <string.h>  // strncpy
#include <stdlib.h>  // rand
#include <math.h> 
#include <limits.h>   // sqrt, pow
//#include <omp.h>     // OpenMP
#include "Timer.h"
#include "Trip.h"


using namespace std;

// Already implemented. see the actual implementations below
void initialize( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] );
void select( Trip trip[CHROMOSOMES], Trip parents[TOP_X] );
void populate( Trip trip[CHROMOSOMES], Trip offsprings[TOP_X] );

// need to implement for your program 1
extern void evaluate( Trip trip[CHROMOSOMES], int coordinates[CITIES][2], float distanceMatrix[CITIES][CITIES] );
extern void crossover( Trip parents[TOP_X], Trip offsprings[TOP_X], float distanceMatrix[CITIES][CITIES], int coordinates[CITIES][2] );
extern void mutate( Trip offsprings[TOP_X] );
extern void setDistanceMatrix(float distanceMatrix[CITIES][CITIES], int coordinates[CITIES][2]);
extern int calIndex(char c);
extern void setFirstOffsprings(float distanceMatrix[CITIES][CITIES], Trip parents[TOP_X], Trip offsprings[TOP_X], int i, int j, bool flag[CITIES]);
extern void setComplement(Trip offsprings[TOP_X], int i);
extern void calTwoDistance(float distanceMatrix[CITIES][CITIES], Trip parents[TOP_X], Trip offsprings[TOP_X], int i, int j, bool flag[CITIES], 
        int first_parent_index,int second_parent_index, int& first_dis, int& second_dis);
extern char randomPick(float distanceMatrix[CITIES][CITIES], Trip offsprings[TOP_X], int i, int j, bool flag[CITIES]);
extern int findIndex(char c, char cArray[CITIES + 1]);
extern bool isValid(char c, bool flag[CITIES]);


/*
 * MAIN: usage: Tsp #threads
 */
int main( int argc, char* argv[] ) {
  Trip trip[CHROMOSOMES];       // all 50000 different trips (or chromosomes)
  Trip shortest;                // the shortest path so far
  int coordinates[CITIES][2];   // (x, y) coordinates of all 36 cities:
  int nThreads = 1;
  float distanceMatrix[CITIES][CITIES];
  
  // verify the arguments
  if ( argc == 2 )
    nThreads = atoi( argv[1] );
  else {
    cout << "usage: Tsp #threads" << endl;
    if ( argc != 1 )
      return -1; // wrong arguments
  }
  cout << "# threads = " << nThreads << endl;

  // shortest path not yet initialized
  shortest.itinerary[CITIES] = 0;  // null path
  shortest.fitness = -1.0;         // invalid distance

  // initialize 50000 trips and 36 cities' coordinates
  initialize( trip, coordinates );

  // start a timer 
  Timer timer;
  timer.start( );

  // change # of threads
 // omp_set_num_threads( nThreads );
  
 // float distanceMatrix[CITIES][CITIES];
  setDistanceMatrix(distanceMatrix, coordinates);
  // find the shortest path in each generation
  for ( int generation = 0; generation < MAX_GENERATION; generation++ ) {

    // evaluate the distance of all 50000 trips
    evaluate(trip, coordinates, distanceMatrix);

    // just print out the progress
    if ( generation % 20 == 0 )
      cout << "generation: " << generation << endl;

    // whenever a shorter path was found, update the shortest path
    if ( shortest.fitness < 0 || shortest.fitness > trip[0].fitness ) {

      strncpy( shortest.itinerary, trip[0].itinerary, CITIES + 1 );
      shortest.fitness = trip[0].fitness;

      cout << "generation: " << generation 
	   << " shortest distance = " << shortest.fitness
	   << "\t itinerary = " << shortest.itinerary << endl;
    }

    // define TOP_X parents and offsprings.
    Trip parents[TOP_X], offsprings[TOP_X];

    // choose TOP_X parents from trip
    select( trip, parents );

    // generates TOP_X offsprings from TOP_X parenets
    crossover( parents, offsprings, distanceMatrix, coordinates );

    // mutate offsprings
    mutate( offsprings );

    // populate the next generation.
    populate( trip, offsprings );
  }

  // stop a timer
  cout << "elapsed time = " << timer.lap( ) << endl;
  return 0;
}
void setDistanceMatrix(float distanceMatrix[CITIES][CITIES], int coordinates[CITIES][2]) {
    int dis = 0;
    for (int i = 0; i < 36; i++) {
        for (int j = 0; j <= i; j++) {
            dis = sqrt(pow(coordinates[i][0] - coordinates[j][0], 2) + pow(coordinates[i][1] - coordinates[j][1], 2));
            distanceMatrix[i][j] = dis;
            distanceMatrix[j][i] = dis;
        }
    }
    return;
}

int compare(const void* a, const void* b) {
    // cout << "entering compare function" << endl;
    float res = ((Trip*) a)->fitness - ((Trip*) b)->fitness;
    if (res > 0) {
        return 1;
    } else if (res < 0) {
        return -1;
    } else {
        return 0;
    }
}

void evaluate(Trip trip[CHROMOSOMES], int coordinates[CITIES][2], float distanceMatrix[CITIES][CITIES]) {
    for (int i = 0; i < CHROMOSOMES; i++) {
        trip[i].fitness = 0;
        for (int j = 0; j < CITIES; j++) {
            if (j == 0) {
                int findex = calIndex(trip[i].itinerary[0]);
                trip[i].fitness += sqrt(pow(coordinates[findex][0], 2) + pow(coordinates[findex][1], 2));
            } else {
                trip[i].fitness += distanceMatrix[calIndex(trip[i].itinerary[j - 1])][calIndex(trip[i].itinerary[j])];
            }
        }
    }
  qsort(trip, CHROMOSOMES, sizeof(Trip), compare);
  
  // for debugging
  /*
  cout << "finishing qsort, the first 10 trip are";
  for (int i = 0; i < 10; i++) {
      cout << "trip[" << i << "] is " << trip[i].itinerary << "," << trip[i].fitness << endl;
  }
  cout << endl;
  */
  return;
}


void crossover(Trip parents[TOP_X], Trip offsprings[TOP_X], float distanceMatrix[CITIES][CITIES], int coordiantes[CITIES][2]) {

    for (int i = 0; i < TOP_X; i += 2) {

       bool flag[CITIES] = {false};
       for (int j = 0; j < CITIES; j++) {
           setFirstOffsprings(distanceMatrix, parents, offsprings, i, j, flag);
       }
       //setFitness(distanceMatrix, offsprings, coordinates, i);
       // for debug, can check if flag has all been set true;
       setComplement(offsprings, i);
    }
    // for debug
    /*
    for (int i = 0; i < 10; i++) {
        cout << "offsprings[" << i << "].itinerary " << offsprings[i].itinerary << "," << offsprings[i].fitness << endl;
    }
    */
}

/*
void setFitness(float distanceMatrix[CITIES][CITIES], Trip offsprings[TOP_X], int coordinates[CITIES][2], int i) {
    int first_city = calIndex(offsprings[i].itinerary[1]);
    offsprings[i].fitness += sqrt(pow(coordinates[first_city][0], 2) + pow(coordinates[first_city][1], 2));
    int second_city;
    for (int i = 2; i <= CITIES; i++) {
        second_city = calIndex(offsprings[i].itinerary[i]);
        offsprings[i].fitness += distanceMatrix[first_city][second_city];
        first_city = second_city;
    }
}*/

void setComplement(Trip offsprings[TOP_X], int index) {
    char benchmark[36] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H','I', 'J', 'K', 'L', 'M', 'N', 'O','P', 'Q',
    'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };
    int orig;
    int com;
    for (int i = 0; i < CITIES; i++) {
        orig = calIndex(offsprings[index].itinerary[i]);
        com = 35 - orig;
        offsprings[index + 1].itinerary[i] = benchmark[com];
    }
    return;
}

void setFirstOffsprings(float distanceMatrix[CITIES][CITIES], Trip parents[TOP_X], Trip offsprings[TOP_X],int i, int j, bool flag[CITIES]) {
    if (j == 0) {
        offsprings[i].itinerary[j] = parents[i].itinerary[j];
        int flag_index = calIndex(parents[i].itinerary[j]);
        flag[flag_index] = true;
        return;
    } 
    int first_dis = 0;
    int second_dis = 0;
    int first_parent_index = findIndex(offsprings[i].itinerary[j - 1], parents[i].itinerary);
    int second_parent_index = findIndex(offsprings[i].itinerary[j - 1], parents[i + 1].itinerary);
    //cout << first_parent_index << "," << second_parent_index << endl;
    calTwoDistance(distanceMatrix, parents, offsprings, i, j, flag, first_parent_index, second_parent_index, first_dis, second_dis );
    
    if (first_dis && second_dis) {
        offsprings[i].itinerary[j] = first_dis < second_dis ? parents[i].itinerary[first_parent_index + 1] : parents[i + 1].itinerary[second_parent_index + 1];
    } else if (first_dis) {
        offsprings[i].itinerary[j] = parents[i].itinerary[first_parent_index + 1];
    } else if (second_dis) {
        offsprings[i].itinerary[j] = parents[i + 1].itinerary[second_parent_index + 1];
    } else {
        offsprings[i].itinerary[j] = randomPick(distanceMatrix, offsprings, i, j, flag);
    }
    flag[calIndex(offsprings[i].itinerary[j])] = true;
    return;
}

void calTwoDistance(float distanceMatrix[CITIES][CITIES], Trip parents[TOP_X], Trip offsprings[TOP_X], int i, int j, bool flag[CITIES],int first_parent_index, 
        int second_parent_index, int& first_dis, int& second_dis) {
   // int first_parent_index = findIndex(offsprings[i].itinerary[j - 1], parents[i].itinerary);
   //int second_parent_index = findIndex(offsprings[i].itinerary[j - 1], parents[i + 1].itinerary);
    int matrix_index_city1 = calIndex(offsprings[i].itinerary[j - 1]);
    if (first_parent_index < CITIES- 1 && isValid(parents[i].itinerary[first_parent_index + 1], flag)) {
        int matrix_index_city2 = calIndex(parents[i].itinerary[first_parent_index + 1]);
        first_dis = distanceMatrix[matrix_index_city1][matrix_index_city2];
    }
    if (second_parent_index < CITIES - 1 && isValid(parents[i + 1].itinerary[second_parent_index + 1], flag)) {
        int matrix_index_city3 = calIndex(parents[i + 1].itinerary[second_parent_index = 1]);
        second_dis = distanceMatrix[matrix_index_city1][matrix_index_city3];
    }
   // cout << "output first and second dis" << first_dis << "and" << second_dis << endl;
    return;
}

bool isValid(char c, bool flag[CITIES]) {
    //cout << "isValid is called" << endl;
    int flag_index = calIndex(c);
    return flag[flag_index] == false;
}

int findIndex(char c, char cArray[CITIES + 1]) {
    //cout << "findIndex is called" << endl;
    char* p = find(cArray, cArray + CITIES, c);
    return p - cArray;
}

int calIndex(char c) {
    //cout << "calIndex is called" << endl;
    return c >= 'A' ? c - 'A' : c - '0' + 26;
}

char randomPick(float distanceMatrix[CITIES][CITIES], Trip offsprings[TOP_X], int i, int j, bool flag[CITIES]) {
    //cout << "randomPick is called" << endl;
    int dis = INT_MAX;
    int city1_index = calIndex(offsprings[i].itinerary[j - 1]);
    int city2_index;
    for (int i = 0; i < CITIES; i++) {
        if (flag[i]) {
            continue;
        }
        if (dis > distanceMatrix[city1_index][i]) {
            dis = distanceMatrix[city1_index][i];
            city2_index = i;
        }
    }
    return city2_index > 25 ? city2_index - 26 + '0' : city2_index + 'A';

}

void mutate(Trip offsprings[TOP_X]) {
    cout << "mutate is called " << endl;
    for (int i = 0; i < TOP_X; i++) {
        if (rand()%100 < 25) {
            int city1 = rand()%36;
            int city2 = rand()%36;
            while (city1 == city2){
                city2 = rand()%36;
            }
            swap(offsprings[i].itinerary[city1], offsprings[i].itinerary[city2]);
        }
    }
}
/*
 * Initializes trip[CHROMOSOMES] with chromosome.txt and coordiantes[CITIES][2] with cities.txt
 *
 * @param trip[CHROMOSOMES]:      50000 different trips
 * @param coordinates[CITIES][2]: (x, y) coordinates of 36 different cities: ABCDEFGHIJKLMNOPQRSTUVWXYZ
 */
void initialize( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] ) {
  // open two files to read chromosomes (i.e., trips)  and cities
  ifstream chromosome_file( "chromosome.txt" );
  ifstream cities_file( "cities.txt" );
  
  // read data from the files
  // chromosome.txt:                                                                                           
  //   T8JHFKM7BO5XWYSQ29IP04DL6NU3ERVA1CZG                                                                    
  //   FWLXU2DRSAQEVYOBCPNI608194ZHJM73GK5T                                                                    
  //   HU93YL0MWAQFIZGNJCRV12TO75BPE84S6KXD
  for ( int i = 0; i < CHROMOSOMES; i++ ) {
    chromosome_file >> trip[i].itinerary;
    trip[i].fitness = 0.0;
  }

  // cities.txt:                                                                                               
  // name    x       y                                                                                         
  // A       83      99                                                                                        
  // B       77      35                                                                                        
  // C       14      64                                                                                        
  for ( int i = 0; i < CITIES; i++ ) {
    char city;
    cities_file >> city;
    int index = ( city >= 'A' ) ? city - 'A' : city - '0' + 26;
    cities_file >> coordinates[index][0] >> coordinates[index][1];
  }

  // close the files.
  chromosome_file.close( );
  cities_file.close( );

  // just for debugging
  if ( DEBUG ) {
    for ( int i = 0; i < CHROMOSOMES; i++ )
      cout << trip[i].itinerary << endl;
    for ( int i = 0; i < CITIES; i++ )
      cout << coordinates[i][0] << "\t" << coordinates[i][1] << endl;
  }
}

/*
 * Select the first TOP_X parents from trip[CHROMOSOMES]
 *
 * @param trip[CHROMOSOMES]: all trips
 * @param parents[TOP_X]:    the firt TOP_X parents
 */
void select( Trip trip[CHROMOSOMES], Trip parents[TOP_X] ) {
  // just copy TOP_X trips to parents
  for ( int i = 0; i < TOP_X; i++ )
    strncpy( parents[i].itinerary, trip[i].itinerary, CITIES + 1 );
}

/*
 * Replace the bottom TOP_X trips with the TOP_X offsprings
 */
void populate( Trip trip[CHROMOSOMES], Trip offsprings[TOP_X] ) {
  // just copy TOP_X offsprings to the bottom TOP_X trips.
  for ( int i = 0; i < TOP_X; i++ )
    strncpy( trip[ CHROMOSOMES - TOP_X + i ].itinerary, offsprings[i].itinerary, CITIES + 1 );

  // for debugging
  if ( DEBUG ) {
    for ( int chrom = 0; chrom < CHROMOSOMES; chrom++ ) 
      cout << "chrom[" << chrom << "] = " << trip[chrom].itinerary 
	   << ", trip distance = " << trip[chrom].fitness << endl;
  }
}
