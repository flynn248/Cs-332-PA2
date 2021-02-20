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
void displayMenu(std::string*, double**, int&);
void doPart(std::string, std::string, std::string);
double** readAdjancyMatrix(std::fstream&, const int);
double* readHeuristicInFile(std::fstream&, const int);
std::string* readNodesInFile(std::fstream &, const int);
std::vector<int> AStar(double[], int&, int&, int**, const int&);

int main()
{
	doPart("PA2 Part 1 Names.csv", "PA2 Part 1 Heuristics.csv", "PA2 Part 1 Distances.csv"); //Part 1
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

	std::fstream readHeuristics(nodeHeuristicFileName, std::ios::in);
	openFileCheck(readHeuristics);

	double* hScore = readHeuristicInFile(readHeuristics, size);

	readHeuristics.close();

	std::fstream readAdjMatrixFile(matrixFileName, std::ios::in);
	openFileCheck(readAdjMatrixFile);

	double** adjMatrix = readAdjancyMatrix(readAdjMatrixFile, size);

	readAdjMatrixFile.close();


	for (int i = 0; i < size; i++)
	{
		std::cout << nodeArray[i] << "   " << hScore[i] << std::endl;
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			std::cout << std::setw(3) << adjMatrix[i][j];
		}
		std::cout << std::endl;
	}


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

void displayMenu(std::string* labelA, double** edgeA, int& size) //display stuff to user
{
	std::cout << "Necessary files have been read!\n" << std::endl;

	std::cout << "\tUniform Cost Search Results :\n"
		<< "Path: " << std::endl
		<< "Cost: " << std::endl
		<< "Explored: " << std::endl
		<< "Fringe: " << std::endl;

	std::cout << "A * Search Results: " << std::endl;
}

std::vector<int> AStar(double hScore[], int& source, int& target, int **adjMatx, const int& SIZE) {
	std::vector<int> openSet = {}; //also called firnge
	std::vector<int> closedSet = {};

	int* prev = new int[SIZE](); //store the previous

	double* gScore = new double[SIZE]; //Cheapest path from source to given node

	double* fScore = new double[SIZE]; //f(n) = g(n) + h(n)

	for (int i = 0; i < SIZE; i++){ //initialize the arrays
		gScore[i] = std::numeric_limits<double>::max();
		fScore[i] = std::numeric_limits<double>::max();
		prev[i] = -1;
	}
	
	gScore[source] = 0;
	fScore[source] = gScore[source] + hScore[source];

	openSet.push_back(source);

	while (openSet.empty() == false){
		
		int currIndex = 0;
		int current = *openSet.begin();

		for (std::vector<int>::iterator it = openSet.begin(); it != openSet.end(); it++, currIndex++) { //current = node in openSet with lowest fScore
			if (fScore[*it] <= fScore[current])
				current = *it;
		}

		openSet.erase(openSet.begin()+currIndex); //remove current from open set
		closedSet.push_back(current); // add to closed set

		if (current == target) {
			std::vector<int> path;
			while (current != -1) {
				path.push_back(current);
				current = prev[current];
			}
			
			std::reverse(path.begin(), path.end());
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
				}

				std::vector<int>::iterator inSet = std::find(closedSet.begin(), closedSet.end(), neighbor);
				if (inSet != closedSet.end()) {
					openSet.erase(inSet);
					openSet.push_back(neighbor);
				}
				else if (!(std::find(openSet.begin(), openSet.end(), neighbor) != openSet.end())) //if neighbor is not in open set
					openSet.push_back(neighbor);
			}
		}
	}

	return;
}