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



/*
Set<int> UniformCostSearch( int source , int target ){
	//A set to store nodes which we have discovered but not explored
	//Sometimes we call this the 'fringe' instead of 'open set'
	Set<int> openSet;

	//A set to store nodes which we have already explored
	Set<int> closedSet;

	//Store the optimal path cost for each node
	//Each index 'i' represents the current minimum path cost for node 'i'
	double[] distances;

	//Store the previous node for each node
	//Each index 'i' represents the current predecessor for node 'i'
	int[] prev;

	//Add the source to the open set
	openSet.Add( source );

	//The shortest path cost from source to source is zero
	distances[ source ] = 0;

	//Keep exploring nodes in the open set and discovering new ones
	while( openSet.Empty() == false ){
		int current = node in openSet with lowest distance;

		//Remove the current node from the open set...
		openSet.Remove( current );

		//...And add it to the closed set so we don't explore it again
		closedSet.Add( current );

		//When we find the target, it's guaranteed to be the shortest path
		//We can safely reconstruct the path and terminate the algorithm
		if( current == target ){
			//We know where each node came from, but not where it goes to
			//So if we want a path from the source to the target...
			//...We have to go from the target to the source then reverse the path

			Set<int> path;
			while( current != -1 ){
				path.Add( current );

				//By setting current to its predecessor, we'll crawl back through the path...
				//...Until we reach the source (which should have a predecessor of -1)
				current = prev[ current ];
			}

			//Since we went from the target to the source, we have to reverse the path
			path.Reverse();
			return path;
		}

		//Discover all the neighbors of current
		foreach( int neighbor of current ){
			//Skip nodes we already explored
			if( closedSet.Contains( neighbor ) )
				continue;

			//Add new nodes to the open set
			if( openSet.Contains( neighbor ) == false )
				openSet.Add( neighbor );

			//Get the edge weight from current to neighbor (from the adjacency matrix)
			double edgeWeight = adj[ current ][ neighbor ];

			//Estimate the cost to go to neighbor from current
			//This is the current minimum path cost to get from source to currrent (distances[ current ])...
			//...Plus the cost to get from current to neighbor (edgeWeight)
			double tentativeDistance = distances[ current ] + edgeWeight;

			//If this tentative cost is better than the one that already exists
			if( tentativeDistance < distances[ neighbor ] ){
				//Update the distance to the new minimum path cost
				distances[ neighbor ] = tentativeDistance;

				//Update the parent node of neighbor
				prev[ neighbor ] = current;
			}
		}
	}

	//If we get to this line then we exhausted all nodes in the open set but never found the target
	//No path exists
	return null;
}
*/


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

/*
Set<int> AStar(double[] hScore, int source, int target) {
	//A set to store nodes which we have discovered but not explored
	//Sometimes we call this the 'fringe' instead of 'open set'
	Set<int> openSet = {};

	//A set to store nodes which we have already explored
	Set<int> closedSet = {};

	//Store the previous node for each node
	//Each index 'i' represents the current predecessor for node 'i'
	int[] prev = {};

	//The 'g' score of each node is the cheapest path from 'source' to the given node
	//This is the same as the 'dist' array in Dijkstra's and UCS
	double[] gScore = {};

	//The 'f' score of each node is g(n) + h(n)
	//This is the evaluation function unique to A*
	//Rather than explore the node on the open set with the lowest gScore (Dijkstra's & UCS)...
	//..A* will explore the node on the open set with the lowest fScore
	//hScore never changes, but whenever gScore[n] is changed, so must fScore[n]
	double[] fScore = {};

	//Initially, we don't know the optimal path cost or the parent of any node
	//To represent a missing node, we'll use -1 as a placeholder
	foreach(int node in graph) {
		gScore[node] = infinity;
		fScore[node] = infinity;
		prev[node] = -1;
	}

	//We know the cost to get to source from source is 0
	//So 0 becomes the gScore of source
	//Since we change gScore, we must also change fScore
	gScore[source] = 0;
	fScore[source] = gScore[source] + hScore[source];

	//Add the source to the open set
	openSet.Add(source);
	//Keep exploring nodes in the open set and discovering new ones
	while (openSet.Empty() == false) {
		int current = node in openSet with lowest fScore;

		//Remove the current node from the open set...
		openSet.Remove(current);

		//...And add it to the closed set so we don't explore it again
		closedSet.Add(current);

		//When we find the target, it's guaranteed to be the shortest path
		//We can safely reconstruct the path and terminate the algorithm
		if (current == target) {
			//We know where each node came from, but not where it goes to
			//So if we want a path from the source to the target...
			//...We have to go from the target to the source then reverse the path

			Set<int> path;
			while (current != -1) {
				path.Add(current);

				//By setting current to its predecessor, we'll crawl back through the path...
				//...Until we reach the source (which should have a predecessor of -1)
				current = prev[current];
			}

			//Since we went from the target to the source, we have to reverse the path
			path.Reverse();
			return path;
		}

		//Discover all the neighbors of current
		foreach(int neighbor of current) {
			//Get the edge weight from current to neighbor (from the adjacency matrix)
			double edgeWeight = adj[current][neighbor];

			//Estimate the cost to go to neighbor from current
			//This is the current minimum path cost to get from source to currrent (gScore[ current ])...
			//...Plus the cost to get from current to neighbor (edgeWeight)
			double tentativeGScore = gScore[current] + edgeWeight;

			//If this tentative cost is better than the one that already exists
			if (tentativeGScore < gScore[neighbor]) {
				//Update the distance to the new minimum path cost
				gScore[neighbor] = tentativeGScore;

				//Since we change gScore, we must also change fScore
				fScore[neighbor] = gScore[neighbor] + hScore[neighbor];

				//Update the parent node of neighbor
				prev[neighbor] = current;

*/
/*
				//If we have found a better path to neighbor, then it must be explored again
				if (closedSet.Contains(neighbor)) {
					closedSet.Remove(neighbor);
					openSet.Add(neighbor);
				}
				//If we are seeing this node for the first time, then ensure that it's on the open set
				else if (openSet.Contains(neighbor) == false)
					openSet.Add(neighbor);
			}
		}
	}

	//If we get to this line then we exhausted all nodes in the open set but never found the target
	//No path exists
	return null;
}
*/