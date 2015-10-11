//
//  SuperGenetic.cpp
//  SuperGenetic
//
//  Created by D, Vivek on 17/06/15.
//  Copyright (c) 2015 Zero. All rights reserved.
//

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <iterator>
#include <map>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <spdlog.h>
#include <vector>

namespace SuperGenetic {
    
	namespace spd = spdlog;
    struct TSPSolver {
        const int CHROMOSOME_SIZE;
        const int POPULATION_SIZE;
        const int MUTATION_POINTS;
        const int MAX_SUBSECTIONS;
		const std::map<int, std::array<int, 2 >> COORDINATES;
		const std::string LOG_FILENAME = std::string("C:\\Users\\vivek\\Documents\\Visual Studio 2015\\Projects\\SuperGenetic\\log");
		std::shared_ptr<spd::logger> logger;
		std::map<std::array<int, 2>, float> distances;

        std::vector<std::vector<int>> population;
		std::vector<int> best_chromosome;
        
        // The constructor that initialitzes all constants
        TSPSolver(const std::map<int, std::array<int, 2 >> coordinates, int chromosome_size,
			int population_size, int mutation_points, int max_subsections) :
			COORDINATES(coordinates), CHROMOSOME_SIZE(chromosome_size), POPULATION_SIZE(population_size),
			MUTATION_POINTS(mutation_points), MAX_SUBSECTIONS(max_subsections) {

			logger = spd::rotating_logger_mt("file_logger", LOG_FILENAME, 1048576 * 5, 3);
			logger->set_pattern("%v");
            logger->info("") <<"City = "<< CHROMOSOME_SIZE<<"\n";
            logger->info("") <<"population = "<< POPULATION_SIZE<<"\n";

			calculate_inter_city_distances();
            this->populate();
            //this->print_population();
        }
        
		void calculate_inter_city_distances() {
			for (int x = 0; x < CHROMOSOME_SIZE - 1; ++x) {
				for (int y = x + 1; y < CHROMOSOME_SIZE; ++y) {
					std::array<int, 2> cityX_coordinates = COORDINATES.find(x + 1)->second;
					std::array<int, 2> cityY_coordinates = COORDINATES.find(y + 1)->second;
					auto cities = std::array<int, 2> {x, y};
					distances[cities] = round(std::sqrt(std::pow(cityX_coordinates[0] - cityY_coordinates[0], 2) +
						std::pow(cityX_coordinates[1] - cityY_coordinates[1], 2)));
					//distances[cities] = psuedo_euclidian_distance(cityX_coordinates[0], cityX_coordinates[1],
					//	cityY_coordinates[0], cityY_coordinates[1]);
				}
			}
		}

		int psuedo_euclidian_distance(const int& x1,
			const int& x2,
			const int& y1,
			const int& y2) {
			int dij = 0;

			int xd = x1 - x2;
			int yd = y1 - y2;

			double rij = sqrt((xd*xd + yd*yd) / 10.0);

			int tij = round(rij);

			if (tij < rij)
			{
				dij = tij + 1;
			}
			else
			{
				dij = tij;
			}

			return dij;
		}

		float distance(int cityX, int cityY) {
			auto city_pair = std::array<int, 2> {cityX, cityY};
			std::sort(city_pair.begin(), city_pair.end());
			return distances[city_pair];
		}

        void populate() {
            std::random_device rd;
            std::mt19937 random_engine(rd());
            
			std::vector<int> base_chromosome(CHROMOSOME_SIZE);
            // fill city_list from from 0 to CHROMOSOME_SIZE - 1
            std::iota(base_chromosome.begin(), base_chromosome.end(), 0);
            
            // Create Chromosomes by permutating the list of cities
            for (auto i = 0; i < POPULATION_SIZE; i++) {
                std::shuffle(base_chromosome.begin(), base_chromosome.end(), random_engine);
                population.push_back(base_chromosome);
            }
            best_chromosome = population.at(0);
        }
        
        void print_population() {
            for(auto chromosome : population) {
                print_chromosome(chromosome);
                logger->info("")<<" F = "<<fitness(chromosome.begin(), chromosome.end())<<"\n";
            }
        }
        
        void print_chromosome(const std::vector<int> chromosome) {
			std::string output;
            for(auto gene : chromosome) {
				output.append(std::to_string(gene));
				output.append("|");
            }
			logger->info(output);
        }
        
        void solve() {
			int i = 0;
            do {
				logger->info("") << "\nGeneration " << i++;
				std::cout << "\nGeneration " << i;
                crossover();
                //logger->info("")<<"\nNew Generation :\n";
                //print_population();
				//logger->info("") << "\nCrossover done. Mutating.....";
                mutate();
				//logger->info("") << "\nMutaion done..";
                supervised_mutate();
				//logger->info("") << "\nSupervised Mutaion done..";
				print_chromosome(best_chromosome);
				const float best_distance = fitness(best_chromosome.begin(), best_chromosome.end());
				logger->info("") << "\n Distance is " << best_distance;
				std::cout << "\n Distance is " << best_distance;
			} while (fitness(best_chromosome.begin(), best_chromosome.end()) > CHROMOSOME_SIZE);
        }
        
        void crossover() {
            std::vector<std::vector<int>> new_generation;
            std::random_device rd;
            std::mt19937 random_engine(rd());
            std::uniform_int_distribution<> population_distribution(0, POPULATION_SIZE - 1);
            std::uniform_int_distribution<> chromosome_distribution(0, CHROMOSOME_SIZE - 1);
            
            for (int i = 0; i < POPULATION_SIZE; i++) {
                // Select 2 chromosomes
				std::vector<int> X = population.at(population_distribution(random_engine));
				std::vector<int> Y = population.at(population_distribution(random_engine));
                
                //logger->info("")<<"\nSelected X = ";print_chromosome(X);
                //logger->info("")<<"\nSelected Y = ";print_chromosome(Y);
                
                // Overlap both chromosome if possible
                int point_of_overlap = chromosome_distribution(random_engine);
                if (chromosomes_can_be_overlapped(X, Y, point_of_overlap) == false) {
                    // Can
                    i--;
                    continue;
                }
                
                //logger->info("")<<"\nOverlapping Chromosomes at position : "<< point_of_overlap <<"\n";
                
				std::vector<int> X_part1(point_of_overlap), X_part2(CHROMOSOME_SIZE - point_of_overlap), Y_part1(point_of_overlap), Y_part2(CHROMOSOME_SIZE - point_of_overlap);
                // Partition X into 2 parts
                std::partition_copy(X.begin(), X.end(), X_part1.begin(), X_part2.begin(),
                                    [=] (int i) {return i < point_of_overlap;});
                // Partition Y into 2 parts
                std::partition_copy(Y.begin(), Y.end(), Y_part1.begin(), Y_part2.begin(),
                                    [=] (int i) {return i < point_of_overlap;});
                
				std::vector<int> child1, child2;
                // Child_1 = X_part1 + Y_part2
                child1.insert(child1.end(), X_part1.begin(), X_part1.end());
                child1.insert(child1.end(), Y_part2.begin(), Y_part2.end());
                
                // Child_2 = Y_part1 + X_part2
                child2.insert(child2.end(), Y_part1.begin(), Y_part1.end());
                child2.insert(child2.end(), X_part2.begin(), X_part2.end());
                
                //logger->info("")<<"\nOverlapping Complete. Selecting Best Chromosome... !\n";
				std::vector<int> current_best_chromosome = find_best_chromosome(std::vector<std::vector<int>>{X,Y,child1,child2});
                // Add best to new generation
                new_generation.push_back(current_best_chromosome);
                
                // Also check & update our population best chromosome
                if (fitness(current_best_chromosome.begin(), current_best_chromosome.end()) <
                    fitness(best_chromosome.begin(), best_chromosome.end())) {
                    best_chromosome = current_best_chromosome;
                }
            }
            
            this->population = new_generation;
            //logger->info("")<<"\nBest Chromosome = ";print_chromosome(best_chromosome);
        }
        
        void mutate() {
            std::random_device rd;
            std::mt19937 random_engine(rd());
            std::uniform_int_distribution<> population_distribution(0, POPULATION_SIZE - 1);
            std::uniform_int_distribution<> chromosome_distribution(0, CHROMOSOME_SIZE - 1);
            
            for(int i = 0; i < MUTATION_POINTS; i++) {
				std::vector<int> X = population.at(population_distribution(random_engine));
                
                // Swap adjacent Genes.
                int gene_to_mutate = chromosome_distribution(random_engine);

                int temp = X[gene_to_mutate];
                X[gene_to_mutate] = X[(gene_to_mutate+1)%CHROMOSOME_SIZE];
                X[(gene_to_mutate+1)%CHROMOSOME_SIZE] = temp;

                // // Swap the 2 randomly selected Genes.
                // int gene_to_mutate1 = chromosome_distribution(random_engine);
                // int gene_to_mutate2 = chromosome_distribution(random_engine);
                
                // int temp = X[gene_to_mutate1];
                // X[gene_to_mutate1] = X[gene_to_mutate2];
                // X[gene_to_mutate2] = temp;

                // Update poulation's best chromosome

                if (fitness(X.begin(), X.end()) <
                    fitness(best_chromosome.begin(), best_chromosome.end())) {
                    best_chromosome = X;
                }
            }
            //logger->info("")<<"\nBest Chromosome = ";print_chromosome(best_chromosome);logger->info("")<<" Fitness = "<<fitness(best_chromosome.begin(), best_chromosome.end());
        }
        
		/**
		==Supervised Mutation==
		Here we find a diseased gene sequence (DGS) by comparing the fitness of subsections
		of the best chromosome in a manner similar to binary search and then optimize it and
		replace the DGS with the optimized one. Only 1/4th of the population is given this replacement,
		provided they have the dgs genes in the order. Another 1/4th of the population is
		given the reverse of the optimized chromosome.
		*/
        void supervised_mutate() {
			auto dgs = find_DGS(best_chromosome);
			auto cured_dgs = nearest_neighbour_solver(dgs);
			replace_dgs_with_cure(dgs, cured_dgs);

			// Do the replace with reverse of cured_dgs
			std::reverse(cured_dgs.begin(), cured_dgs.end());
			replace_dgs_with_cure(dgs, cured_dgs);
        }

		/**
		Try and Replace 1/4th of population with cured_dgs
		*/
		void replace_dgs_with_cure(const std::vector<int> dgs, const std::vector<int> cured_dgs) {
			auto dgs_size = dgs.size();
			std::random_device rd;
			std::mt19937 random_engine(rd());
			std::uniform_int_distribution<> population_distribution(0, POPULATION_SIZE - 1);

			for (int i = 0; i < POPULATION_SIZE / 4; ++i) {
				std::vector<int> X = population.at(population_distribution(random_engine));
				auto start = X.begin();
				auto end = X.end();
				std::vector<int> positions;
				for (int j = 0; j < dgs_size; ++j) {
					start = std::find(start, end, dgs[j]);
					positions.push_back(start - X.begin());
					if (start == end) {
						break;
					}
				}
				//if the gene from dgs were not found in the same order, move to another chromosome
				if (start == end) {
					continue;
				}

				for (int j = 0; j < dgs_size; ++j) {
					X[positions[j]] = dgs[j];
				}
			}
		}
        
        std::vector<int> find_DGS(std::vector<int> X) {
            auto start = X.begin();
            auto end = X.end();
            // find Diseased Gene Sequence
            for (int i = 0; i < MAX_SUBSECTIONS; i++) {
                long mid = (end - start)/2;
                float part1_fitness = fitness(start, start + mid);
                float part2_fitness = fitness(start + mid, end);
                if (part1_fitness > part2_fitness) {
                    end = start + mid;
                }
                else {
                    start = start + mid;
                }
            }
			std::vector<int> dgs;
			std::copy(start, end, std::back_inserter(dgs));
            return dgs;
        }
    
		std::vector<int> nearest_neighbour_solver(const std::vector<int> cities) {
			std::vector<int> visited_cities;
			std::set<int> unvisited_cities(cities.begin(), cities.end());
			int current_city = *unvisited_cities.begin();
			visited_cities.push_back(current_city);
			unvisited_cities.erase(current_city);
			for (unsigned int i = 0; i < cities.size() - 1; ++i) {
				current_city = find_nearest_city(current_city, unvisited_cities);
				visited_cities.push_back(current_city);
				unvisited_cities.erase(current_city);
			}
			return visited_cities;
		}

		int find_nearest_city(const int city, const std::set<int> unvisited_cities) {
			assert(unvisited_cities.size() != 0);
			float shortest_distance = 0;
			int nearest_city = *unvisited_cities.begin();
			for (int current_city : unvisited_cities) {
				float current_distance = distance(current_city, city);
				if (current_distance < shortest_distance) {
					shortest_distance = current_distance;
					nearest_city = current_city;
				}
			}
			return nearest_city;
		}

    private:
        bool chromosomes_can_be_overlapped(const std::vector<int> X, const std::vector<int> Y, int point_of_overlap) {
            std::vector<int> intersection1;
            std::vector<int> intersection2;
			std::vector<int> orderedX;
			std::vector<int> orderedY;
			std::copy(X.begin(), X.end(), std::back_inserter(orderedX));
			std::copy(Y.begin(), Y.end(), std::back_inserter(orderedY));
			std::sort(orderedX.begin(), orderedX.end());
			std::sort(orderedY.begin(), orderedY.end());

            std::set_intersection(orderedX.begin(), orderedX.begin() + point_of_overlap, orderedY.begin() + point_of_overlap, orderedY.end(), std::back_inserter(intersection1));
            std::set_intersection(orderedY.begin(), orderedY.begin() + point_of_overlap, orderedX.begin() + point_of_overlap, orderedX.end(), std::back_inserter(intersection2));
            if (intersection1.size() > 0 || intersection2.size() > 0) {
                return false;
            }
            return true;
        }
        
        // Lower the Better
        float fitness(std::vector<int>::iterator start, std::vector<int>::iterator end) {
            float total_distance = 0;
            int first = 0;
            if(start != end) {
                first = *start;
            }
            start++;
            while (start != end) {
				total_distance += distance( *start, first);
                first = *start;
                start++;
            }
            return total_distance;
        }
        
		std::vector<int> find_best_chromosome(const std::vector<std::vector<int>> chromosomes) {
            if (chromosomes.size() == 0) {
                throw std::invalid_argument("Empty list of Chromosomes is not acceptable");
            }
			std::vector<int> best_chromosome = chromosomes.at(0);
            float best_fitness = fitness(best_chromosome.begin(), best_chromosome.end());
            for (auto chromosome : chromosomes) {
                float current_fitness = fitness(chromosome.begin(), chromosome.end());
                // Lower the Better
                if (best_fitness > current_fitness) {
                    best_fitness = current_fitness;
                    best_chromosome = chromosome;
                }
            }
            return best_chromosome;
        }


		
    };

}
