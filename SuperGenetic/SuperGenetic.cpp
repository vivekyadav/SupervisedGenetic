//
//  SuperGenetic.cpp
//  SuperGenetic
//
//  Created by D, Vivek on 17/06/15.
//  Copyright (c) 2015 Zero. All rights reserved.
//

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <vector>

namespace SuperGenetic {
    
    typedef int Gene;
    typedef std::vector<Gene> Chromosome;
    
    struct TSPSolver {
        const int CHROMOSOME_SIZE;
        const int DGS_SIZE;
        const int POPULATION_SIZE;
        const int MUTATION_POINTS;
        const int MAX_SUBSECTIONS;
        std::vector<Chromosome> population;
        const std::map<std::set<int>, float> distances;
        Chromosome best_chromosome;
        
        // The constructor that initialitzes all constants
        TSPSolver(int chromosome_size, int dgs_size, int population_size, int mutation_points, int max_subsections) :
        CHROMOSOME_SIZE(chromosome_size), DGS_SIZE(dgs_size), POPULATION_SIZE(population_size), MUTATION_POINTS(mutation_points), MAX_SUBSECTIONS(max_subsections) {
            std::cout<<"City = "<< CHROMOSOME_SIZE<<"\n";
            std::cout<<"DGS = "<< DGS_SIZE<<"\n";
            std::cout<<"population = "<< POPULATION_SIZE<<"\n";
            this->populate();
            this->print_population();
        }
        
        void populate() {
            std::random_device rd;
            std::mt19937 random_engine(rd());
            
            Chromosome base_chromosome(CHROMOSOME_SIZE);
            // fill city_list from from 0 to CHROMOSOME_SIZE - 1
            std::iota(base_chromosome.begin(), base_chromosome.end(), 0);
            
            // Create Chromosomes by permutating the list of cities
            for (auto i = 0; i < POPULATION_SIZE; i++) {
                std::shuffle(base_chromosome.begin(), base_chromosome.end(), random_engine);
                population.push_back(Chromosome(base_chromosome));
            }
            best_chromosome = population.at(0);
            
            //Calculate distance
            calculate_distances();
        }
        
        void print_population() {
            for(auto chromosome : population) {
                print_chromosome(chromosome);
                std::cout<<" F = "<<fitness(chromosome.begin(), chromosome.end())<<"\n";
            }
        }
        
        void print_chromosome(const Chromosome chromosome) {
            for(auto gene : chromosome) {
                std::cout<<gene<<"|";
            }
        }
        
        void solve() {
            for (int i = 0; i < 5; i++) {
                crossover();
                //std::cout<<"\nNew Generation :\n";
                //print_population();
                mutate();
                //supervised_mutate();
            }
        }
        
        void crossover() {
            std::vector<Chromosome> new_generation;
            std::random_device rd;
            std::mt19937 random_engine(rd());
            std::uniform_int_distribution<> population_distribution(0, POPULATION_SIZE - 1);
            std::uniform_int_distribution<> chromosome_distribution(0, CHROMOSOME_SIZE - 1);
            
            for (int i = 0; i < POPULATION_SIZE; i++) {
                // Select 2 chromosomes
                Chromosome X = population.at(population_distribution(random_engine));
                Chromosome Y = population.at(population_distribution(random_engine));
                
                //std::cout<<"\nSelected X = ";print_chromosome(X);
                //std::cout<<"\nSelected Y = ";print_chromosome(Y);
                
                // Overlap both chromosome if possible
                int point_of_overlap = chromosome_distribution(random_engine);
                if (chromosomes_can_be_overlapped(X, Y, point_of_overlap) == false) {
                    // Can
                    i--;
                    continue;
                }
                
                //std::cout<<"\nOverlapping Chromosomes at position : "<< point_of_overlap <<"\n";
                
                Chromosome X_part1(point_of_overlap), X_part2(CHROMOSOME_SIZE - point_of_overlap), Y_part1(point_of_overlap), Y_part2(CHROMOSOME_SIZE - point_of_overlap);
                // Partition X into 2 parts
                std::partition_copy(X.begin(), X.end(), X_part1.begin(), X_part2.begin(),
                                    [=] (int i) {return i < point_of_overlap;});
                // Partition Y into 2 parts
                std::partition_copy(Y.begin(), Y.end(), Y_part1.begin(), Y_part2.begin(),
                                    [=] (int i) {return i < point_of_overlap;});
                
                Chromosome child1, child2;
                // Child_1 = X_part1 + Y_part2
                child1.insert(child1.end(), X_part1.begin(), X_part1.end());
                child1.insert(child1.end(), Y_part2.begin(), Y_part2.end());
                
                // Child_2 = Y_part1 + X_part2
                child2.insert(child2.end(), Y_part1.begin(), Y_part1.end());
                child2.insert(child2.end(), X_part2.begin(), X_part2.end());
                
                //std::cout<<"\nOverlapping Complete. Selecting Best Chromosome... !\n";
                Chromosome current_best_chromosome = find_best_chromosome(std::vector<Chromosome>{X,Y,child1,child2});
                // Add best to new generation
                new_generation.push_back(current_best_chromosome);
                
                // Also check & update our population best chromosome
                if (fitness(current_best_chromosome.begin(), current_best_chromosome.end()) <
                    fitness(best_chromosome.begin(), best_chromosome.end())) {
                    best_chromosome = current_best_chromosome;
                }
            }
            
            this->population = new_generation;
            //std::cout<<"\nBest Chromosome = ";print_chromosome(best_chromosome);
        }
        
        void mutate() {
            std::random_device rd;
            std::mt19937 random_engine(rd());
            std::uniform_int_distribution<> population_distribution(0, POPULATION_SIZE - 1);
            std::uniform_int_distribution<> chromosome_distribution(0, CHROMOSOME_SIZE - 1);
            
            for(int i = 0; i < MUTATION_POINTS; i++) {
                Chromosome X = population.at(population_distribution(random_engine));
                int gene_to_mutate = chromosome_distribution(random_engine);
                
                // Swap adjacent Genes.
                Gene temp = X[gene_to_mutate];
                X[gene_to_mutate] = X[(gene_to_mutate+1)%CHROMOSOME_SIZE];
                X[(gene_to_mutate+1)%CHROMOSOME_SIZE] = temp;
                
                // Update poulation's best chromosome
                if (fitness(X.begin(), X.end()) <
                    fitness(best_chromosome.begin(), best_chromosome.end())) {
                    best_chromosome = X;
                }
            }
            std::cout<<"\nBest Chromosome = ";print_chromosome(best_chromosome);std::cout<<" Fitness = "<<fitness(best_chromosome.begin(), best_chromosome.end());
        }
        
        void supervised_mutation() {
            auto dgs = find_DGS(best_chromosome);
        }
        
        std::tuple<Chromosome::iterator, Chromosome::iterator> find_DGS(Chromosome X) {
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
            return std::make_tuple(start, end);
        }
        
    private:
        bool chromosomes_can_be_overlapped(Chromosome X, Chromosome Y, int point_of_overlap) {
            std::vector<Gene> intersection1;
            std::vector<Gene> intersection2;
            std::set_intersection(X.begin(), X.begin() + point_of_overlap, Y.begin() + point_of_overlap, Y.end(), std::back_inserter(intersection1));
            std::set_intersection(Y.begin(), Y.begin() + point_of_overlap, X.begin() + point_of_overlap, X.end(), std::back_inserter(intersection2));
            if (intersection1.size() > 0 || intersection2.size() > 0) {
                return false;
            }
            return true;
        }
        
        // Lower the Better
        float fitness(Chromosome::iterator start, Chromosome::iterator end) {
            float total_distance = 0;
            Gene first = 0;
            if(start != end) {
                first = *start;
            }
            start++;
            while (start != end) {
                total_distance += abs(*start - first);
                first = *start;
                start++;
            }
            return total_distance;
        }
        
        Chromosome find_best_chromosome(const std::vector<Chromosome> chromosomes) {
            if (chromosomes.size() == 0) {
                throw std::invalid_argument("Empty list of Chromosomes is not acceptable");
            }
            Chromosome best_chromosome = chromosomes.at(0);
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
        
        void calculate_distances() {
            //assign random postitions to Geness/cities
            //std::vector<std::tuplepositions
        }
    };
}
