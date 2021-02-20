//Developed by Shane Flynn
//Completed using MS Visual Studio 2019
#include <algorithm> //might need to get rid of
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <limits>


void openFileCheck(std::fstream&);
int countNodesInFile(std::fstream&);
void doPart(std::string, std::string, std::string);
double** readAdjancyMatrix(std::fstream&, const int);
double* readHeuristicInFile(std::fstream&, const int);
std::string* readNodesInFile(std::fstream &, const int);
std::vector<int> AStar(double[], int, int, double**, const int&, std::string *);

int main()
{
	std::cout << "Graph #1\n\n";
	doPart("PA2 Part 1 Names.csv", "PA2 Part 1 Heuristics.csv", "PA2 Part 1 Distances.csv"); //Part 1
	std::cout << "Graph #2\n\n";
	doPart("PA2 Part 2 Names.csv", "PA2 Part 2 Heuristics.csv", "PA2 Part 2 Distances.csv"); //Part 2

	std::cin.get();
	return 0;
}

void doPart(std::string nodeNamesFileName, std::string nodeHeuristicFileName, std::string matrixFileName) {
	std::fstream readNames(nodeNamesFileName, std::ios::in);
	openFileCheck(readNames);

	int size = countNodesInFile(readNames);
	readNames.seekg(0); //move cursor to beginning of file
	std::string* nodeArray = readNodesInFile(readNames, size);

	readNames.close();

	std::fstream readAdjMatrixFile(matrixFileName, std::ios::in);
	openFileCheck(readAdjMatrixFile);

	double** adjMatrix = readAdjancyMatrix(readAdjMatrixFile, size);

	readAdjMatrixFile.close();

	double* hScore = new double[size]();

	std::cout << "Uniform Cost Search Results:\n";
	AStar(hScore, 0, size - 1, adjMatrix, size, nodeArray);

	std::fstream readHeuristics(nodeHeuristicFileName, std::ios::in);
	openFileCheck(readHeuristics);
	hScore = readHeuristicInFile(readHeuristics, size);

	readHeuristics.close();

	std::cout << "A* Search Results:\n";
	AStar(hScore, 0, size - 1, adjMatrix, size, nodeArray);


	for (int i = 0; i < size; i++) //deleting row of pointers
		delete[] adjMatrix[i];

	delete[] adjMatrix;
	adjMatrix = nullptr;
	delete[] nodeArray;
	nodeArray = nullptr;
	delete[] hScore;
	hScore = nullptr;
}

int countNodesInFile(std::fstream& file) {

	int size = 0;
	std::string temp;

	while (std::getline(file, temp, ',')) //find total # of nodes in graph
	{
		size++;
		if (file.eof()) break;
	}
	return size;
}

std::string* readNodesInFile(std::fstream& fileNodes, const int size) {

	std::string *nodeArray = new std::string[size];

	for (int i = 0; i < size; i++) //Parse the location values from file into labelArray
	{
		if (i != size - 1)
			std::getline(fileNodes, nodeArray[i], ',');
		else
			std::getline(fileNodes, nodeArray[i], '\n');
	}

	return nodeArray;
}

double* readHeuristicInFile(std::fstream& fileNodes, const int size) {
	
	std::string temp;
	double* nodeArray = new double[size];

	for (int i = 0; i < size; i++) //Parse the location values from file into labelArray
	{
		if (i != size - 1) {
			std::getline(fileNodes, temp, ',');
			nodeArray[i] = std::stod(temp);
		}
		else {
			std::getline(fileNodes, temp, '\n');
			nodeArray[i] = std::stod(temp);
		}
	}

	return nodeArray;
}

double** readAdjancyMatrix(std::fstream& file, const int size) {

	double** edgeArray = new double* [size]; //array to hold edge weights
	std::string temp;

	for (int i = 0; i < size; i++) //making array 2D
	{
		edgeArray[i] = new double[size];
	}

	for (int i = 0; i < size; i++) //parse distances into array
	{
		for (int j = 0; j < size; j++)
		{
			if (j != size - 1) //Prevents a bug from reading the '\n' as an input
			{
				std::getline(file, temp, ',');
				edgeArray[i][j] = std::stod(temp);
			}
			else
			{
				std::getline(file, temp, '\n');
				edgeArray[i][j] = std::stod(temp);
			}
		}
	}

	return edgeArray;

}

void openFileCheck(std::fstream& file) //check if file failed to open
{
	if (!file)
	{
		std::cout << "File failed to open..." << std::endl;
		std::cout << "Terminating program..." << std::endl;
		std::cout << "Hit enter to close window...   " << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
}

std::vector<int> AStar(double hScore[], int source, int target, double **adjMatx, const int& SIZE, std::string *nodeNames) {
	std::vector<int> openSet = {}; //also called firnge
	std::vector<int> closedSet = {};

	int* prev = new int[SIZE](); //store the previous

	double* gScore = new double[SIZE](); //Cheapest path from source to given node

	double* fScore = new double[SIZE](); //f(n) = g(n) + h(n)

	for (int i = 0; i < SIZE; i++){ //initialize the arrays
		gScore[i] = std::numeric_limits<double>::max();
		fScore[i] = std::numeric_limits<double>::max();
		prev[i] = -1;
	}

	gScore[source] = 0;
	fScore[source] = gScore[source] + hScore[source];

	openSet.push_back(source);

	while (openSet.empty() == false){

		int current = *openSet.begin();

		std::vector<int>::iterator curr = openSet.begin();
		for (std::vector<int>::iterator it = openSet.begin(); it != openSet.end(); it++) { //current = node in openSet with lowest fScore
			if (fScore[*it] <= fScore[current]) {
				current = *it;
				curr = it;
			}
		}

		openSet.erase(curr); //remove current from open set
		closedSet.push_back(current); // add to closed set

		if (current == target) {
			std::vector<int> path;

			int cost = gScore[current];

			while (current != -1) {
				path.push_back(current);
				current = prev[current];
			}
			
			std::reverse(path.begin(), path.end());

			//Output of results
			std::cout << "Path: [";
			for (std::vector<int>::iterator it = path.begin(); it != path.end(); it++) {
				if (it == path.end()-1)
					std::cout << nodeNames[*it];
				else
					std::cout << nodeNames[*it] << "->";
			}
			std::cout << "]\n"
				<< "Cost: " << cost << std::endl
				<< "Explored: ";

			for (std::vector<int>::iterator it = closedSet.begin(); it != closedSet.end(); it++) {
				std::cout << nodeNames[*it] << " ";
			}
			std::cout << "\nFringe: ";

			for (std::vector<int>::iterator it = openSet.begin(); it != openSet.end(); it++) {
				std::cout << nodeNames[*it] << " ";
			}
			std::cout << std::endl << std::endl;


			return path;
		}

		for (int neighbor = 0; neighbor < SIZE; neighbor++) {
			if (adjMatx[current][neighbor] != -1) {
				double edgeWeight = adjMatx[current][neighbor];

				double tentativeGScore = gScore[current] + edgeWeight;

				if (tentativeGScore < gScore[neighbor]) {
					gScore[neighbor] = tentativeGScore;

					fScore[neighbor] = gScore[neighbor] + hScore[neighbor];

					prev[neighbor] = current;


					std::vector<int>::iterator inSet = std::find(closedSet.begin(), closedSet.end(), neighbor);
					if (inSet != closedSet.end()) {
						closedSet.erase(inSet);
						openSet.push_back(neighbor);
					}
					else if (!(std::find(openSet.begin(), openSet.end(), neighbor) != openSet.end())) //if neighbor is not in open set
						openSet.push_back(neighbor);
				}
			}
		}
	}
	
	return std::vector<int>(); //Return null if no path found
}