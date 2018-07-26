#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/graph/graphviz.hpp>
#include <libs/graph/src/read_graphviz_new.cpp>
#include <string>
#include <stack>
#include <set>


/**
 * Vertex Data used to store the information used by the Tarjan version of the algorithm
 * lowpt, lowvine and number are the same as used in the pseudo-code of the algorithm.
 * number is initialized to -1 as a convention to say that a vertex is not yet numbered
 * component is an integer indicating to which of the SCC of the grap the node belongs to after the algorithm is done running.
 */
struct TarjanVertexData {
	int lowpt;
	int lowvine;
	int number = -1;
	int component;
};

/**
* Vertex Data used to store the information used by the Nuutila version of the algorithm
* root points to the root vertex of any given vertex.
* inComponent is the same as in the pseudo-code of the algorithm
* visited is a boolean to keep track of when a node has been visited by the algorithm
* visitIndex is an integer keeping track of the order vertices are visited
*/
struct NuutilaVertexData {
	boost::adjacency_list<>::vertex_descriptor root;
	bool inComponent;
	bool visited = false;
	int visitIndex;
};

/**
* Vertex Data used to store the information used by the Pearce version of the algorithm
* both rindex and root are the same as in the pseudo-code of the algorithm
*/
struct PearceVertexData {
	int rindex;
	bool root = false;
};



/**
* Simple typedef of a Directed Graph
*/
typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::no_property> DirectedGraph;

/**
* typedef of a Directed Graph containing in it's vertices the TarjanVertexData
*/
typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, TarjanVertexData, boost::no_property> TarjanDirectedGraph;

/**
* typedef of a Directed Graph containing in it's vertices the NuutilaVertexData
*/
typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, NuutilaVertexData, boost::no_property> NuutilaDirectedGraph;

/**
* typedef of a Directed Graph containing in it's vertices the PearceVertexData
*/
typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, PearceVertexData, boost::no_property> PearceDirectedGraph;



/**
* Function identifier of createRandomDirectedGraph. Further documentation on the function code body.
*/
DirectedGraph* createRandomDirectedGraph(int numVertex);

/**
* Function identifier of createRandomTarjanDirectedGraph. Further documentation on the function code body.
*/
TarjanDirectedGraph* createRandomTarjanDirectedGraph(int numVertex);

/**
* Function identifier of createRandomNuutilaDirectedGraph. Further documentation on the function code body.
*/
NuutilaDirectedGraph* createRandomNuutilaDirectedGraph(int numVertex);

/**
* Function identifier of createRandomPearceDirectedGraph. Further documentation on the function code body.
*/
PearceDirectedGraph* createRandomPearceDirectedGraph(int numVertex);

/**
* Function identifier of saveDirectedGraphToFile. Further documentation on the function code body.
*/
void saveDirectedGraphToFile(DirectedGraph* g,int numVertex, int index);

/**
* Function identifier of saveDirectedGraphToFile. Further documentation on the function code body.
*/
DirectedGraph* readDirectedGraphFromFile(int numVertex, int index);

/**
* Function identifier of generateAndSaveRandomGraphs. Further documentation on the function code body.
*/
void generateAndSaveRandomGraphs(int numGraphs);



/**
* Function identifier of tarjanSCC. Further documentation on the function code body.
*/
void tarjanSCC(TarjanDirectedGraph* g);

/**
* Function identifier of tarjanStrongConnect. Further documentation on the function code body.
*/
void tarjanStrongConnect(boost::adjacency_list<>::vertex_descriptor* v, 
	TarjanDirectedGraph& g);

/**
* Function identifier of displayTarjanSCC. Further documentation on the function code body.
*/
TarjanDirectedGraph displayTarjanSCC(TarjanDirectedGraph* g);



/**
* Function identifier of nuutilaSCC. Further documentation on the function code body.
*/
void nuutilaSCC(NuutilaDirectedGraph* g);

/**
* Function identifier of nuutilaVisit. Further documentation on the function code body.
*/
void nuutilaVisit(NuutilaDirectedGraph::vertex_descriptor* v,
	NuutilaDirectedGraph& g);


/**
* Function identifier of imperativePearceSCC. Further documentation on the function code body.
*/
void imperativePearceSCC(PearceDirectedGraph* g);

/**
* Function identifier of pearceVisit. Further documentation on the function code body.
*/
void pearceVisit(PearceDirectedGraph* g, 
	PearceDirectedGraph::vertex_descriptor* v, 
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, 
	std::stack<int>* iteratorStack, 
	int* c, 
	int* index);

/**
* Function identifier of pearceBeginVisiting. Further documentation on the function code body.
*/
void pearceBeginVisiting(PearceDirectedGraph* g,
	PearceDirectedGraph::vertex_descriptor* v, 
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, 
	std::stack<int>* iteratorStack, 
	int* c, 
	int* index);

/**
* Function identifier of pearceVisitLoop. Further documentation on the function code body.
*/
void pearceVisitLoop(PearceDirectedGraph* g,
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack,
	std::stack<int>* iteratorStack,
	int* c, 
	int* index);

/**
* Function identifier of pearceFinishVisiting. Further documentation on the function code body.
*/
void pearceFinishVisiting(PearceDirectedGraph* g,
	PearceDirectedGraph::vertex_descriptor* v,
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack,
	std::stack<int>* iteratorStack,
	int* c,
	int* index);

/**
* Function identifier of pearceBeginEdge. Further documentation on the function code body.
*/
bool pearceBeginEdge(PearceDirectedGraph* g,
	PearceDirectedGraph::vertex_descriptor* v,
	int* i, 
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack,
	std::stack<int>* iteratorStack, 
	int* c, 
	int* index);

/**
* Function identifier of pearceFinishEdge. Further documentation on the function code body.
*/
void pearceFinishEdge(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, int* i);


int nuutilaCounter;
std::stack<NuutilaDirectedGraph::vertex_descriptor> nuutilaVerticesStack; //stack as in the pseudo-code
std::set<NuutilaDirectedGraph::vertex_descriptor> nuutilaVerticesSet; // auxiliary set, used for checking if an element is inside the stack

std::stack<TarjanDirectedGraph::vertex_descriptor> tarjanPointsStack; //same stack used in the pesudo-code
std::set<TarjanDirectedGraph::vertex_descriptor> tarjanPointsSet; //auxiliary set, used for checking if an element is inside the stack
int tarjanI; // same i used in the pseudo-code
int tarjanComponentsCount; // counts how many Strongly Connected Components there are in the Graph, used for the generation of the final graph.


/**
 * Simple main function that only calls the creation of a graph, saving of a graph to file, and one of the algorithms.
 */
int main() {
	std::cout << "Advanced Algorithms Project" << std::endl;

	//DirectedGraph* g;
	//g = readDirectedGraphFromFile(5, 0);
	//delete g;
	//PearceDirectedGraph* g = createRandomPearceDirectedGraph(5);
	//imperativePearceSCC(g);
	//NuutilaDirectedGraph *g = createRandomNuutilaDirectedGraph(5);
	//nuutilaSCC(g);
	TarjanDirectedGraph* g = createRandomTarjanDirectedGraph(5);
	TarjanDirectedGraph sccGraph = displayTarjanSCC(g);
	//generateAndSaveRandomGraphs(100);
	std::cin.get();
	return 0;
}

/**
* Creates a random Directed Graph with exactly numVertex vertices.
* Uses an uniform distribution in the range 1 to 10000.
* First, a treshold is set. Then for every possible pair of vertices (including self loops) a new random value is calculated,
* and if the random value is greater than the treshold, the edge is added to the graph.
*/
DirectedGraph* createRandomDirectedGraph(int numVertex) {

	boost::random::mt19937 gen(time(NULL));
	boost::random::uniform_int_distribution<> dist(1, 10000);
	int treshold = dist(gen); //sets treshold
	
	auto g = new DirectedGraph;
	auto vertices = new boost::adjacency_list<>::vertex_descriptor[numVertex];
	for (int i = 0; i < numVertex; i++) {
		vertices[i] = add_vertex(*g);//adds numVertex vertices to the graph
	}
	for (int i = 0, randValue = 0; i < numVertex; i++) {
		for (int j = 0; j < numVertex; j++) {
			randValue = dist(gen); //random value for the edge
			if (randValue >= treshold) { 
				auto e = boost::add_edge(i, j, *g); //adds edge if the value is greater than the treshold
			}
		}
	}
	delete[] vertices;
	return g;
}

/**
* Creates a random Tarjan Directed Graph with exactly numVertex vertices.
* Uses an uniform distribution in the range 1 to 10000.
* First, a treshold is set. Then for every possible pair of vertices (including self loops) a new random value is calculated,
* and if the random value is greater than the treshold, the edge is added to the graph.
*/
TarjanDirectedGraph* createRandomTarjanDirectedGraph(int numVertex) {

	boost::random::mt19937 gen(time(NULL));
	boost::random::uniform_int_distribution<> dist(1, 10000);
	int treshold = dist(gen);

	auto g = new TarjanDirectedGraph;
	auto vertices = new boost::adjacency_list<>::vertex_descriptor[numVertex];
	for (int i = 0; i < numVertex; i++) {
		vertices[i] = add_vertex(*g);
	}
	for (int i = 0, randValue = 0; i < numVertex; i++) {
		for (int j = 0; j < numVertex; j++) {
			randValue = dist(gen);
			if (randValue >= treshold) {
				auto e = boost::add_edge(i, j, *g);
			}
		}
	}

	delete[] vertices;
	return g;

}

/**
* Creates a random Nuutila DIrected Graph with exactly numVertex vertices.
* Uses an uniform distribution in the range 1 to 10000.
* First, a treshold is set. Then for every possible pair of vertices (including self loops) a new random value is calculated,
* and if the random value is greater than the treshold, the edge is added to the graph.
*/
NuutilaDirectedGraph* createRandomNuutilaDirectedGraph(int numVertex) {

	boost::random::mt19937 gen(time(NULL));
	boost::random::uniform_int_distribution<> dist(1, 10000);
	int treshold = dist(gen);

	auto g = new NuutilaDirectedGraph;
	auto vertices = new boost::adjacency_list<>::vertex_descriptor[numVertex];
	for (int i = 0; i < numVertex; i++) {
		vertices[i] = add_vertex(*g);
	}
	for (int i = 0, randValue = 0; i < numVertex; i++) {
		for (int j = 0; j < numVertex; j++) {
			randValue = dist(gen);
			if (randValue >= treshold) {
				auto e = boost::add_edge(i, j, *g);
			}
		}
	}

	delete[] vertices;
	return g;

}

/**
* Creates a random Pearce Directed Graph with exactly numVertex vertices.
* Uses an uniform distribution in the range 1 to 10000.
* First, a treshold is set. Then for every possible pair of vertices (including self loops) a new random value is calculated,
* and if the random value is greater than the treshold, the edge is added to the graph.
*/ 
PearceDirectedGraph* createRandomPearceDirectedGraph(int numVertex) {

	boost::random::mt19937 gen(1);
	boost::random::uniform_int_distribution<> dist(1, 10000);
	int treshold = dist(gen);

	auto g = new PearceDirectedGraph;
	auto vertices = new boost::adjacency_list<>::vertex_descriptor[numVertex];
	for (int i = 0; i < numVertex; i++) {
		vertices[i] = add_vertex(*g);
	}
	for (int i = 0, randValue = 0; i < numVertex; i++) {
		for (int j = 0; j < numVertex; j++) {
			randValue = dist(gen);
			if (randValue >= treshold) {
				auto e = boost::add_edge(i, j, *g);
			}
		}
	}

	delete[] vertices;
	return g;

}

/**
* Saves a Directed Graph into a file, naming it DirectedGraph, followed by the number of vertices, and then an id.
*/
void saveDirectedGraphToFile(DirectedGraph* g, int numVertex, int index) {

	std::string fileName = "DirectedGraph" + std::to_string(numVertex) + "Vertexes" + std::to_string(index);
	std::string path = "./DirectedGraphs/" + fileName + ".txt";
	std::ofstream writer(path);
	writer.close();
	return;
}

DirectedGraph* readDirectedGraphFromFile(int numVertex, int index) {
	std::string fileName = "DirectedGraph" + std::to_string(numVertex) + "Vertexes" + std::to_string(index);
	std::string path = "./DirectedGraphs/" + fileName + ".txt";
	std::ifstream reader(path);
	DirectedGraph* g = new DirectedGraph;
	boost::dynamic_properties dp;

	if (reader) {
		if (boost::read_graphviz(path, *g, dp)) {
			return g;
		}
		else {
			std::cout << "Problem Reading File" << std::endl;
			return g;
		}
	}
	else {
		std::cout << "Problem Reading File" << std::endl;
		return NULL;
	}
}

/**
* Generates numGraphs random graphs and saves all of them to file.
* A fourth of numGraph graphs are generated with 5, 10, 20 and 50 vertexes
*/
void generateAndSaveRandomGraphs(int numGraphs) {
	int numVertex = 5;
	DirectedGraph* g;
	for (int i = 0; i < 4; i++) {
		switch (i) {
		case 1: numVertex = 10;
			break;
		case 2: numVertex = 20;
			break;
		case 3: numVertex = 50;
			break;
		}
		for (int j = 0; j < numGraphs / 4; j++) {
			g = createRandomDirectedGraph(numVertex);
			saveDirectedGraphToFile(g, numVertex, j);
			delete g;
		}
	}

}


/**
* Builds a new Graph based on the result of the Tarjan SCC algorithm result.
*/
TarjanDirectedGraph displayTarjanSCC(TarjanDirectedGraph* g) {
	tarjanSCC(g);
	TarjanDirectedGraph resultGraph(num_vertices(*g));

	auto iterations = boost::edges(*g);
	auto It = iterations.first;
	for (; It != iterations.second; It++) {
		if ((*g)[source(*It, *g)].component == (*g)[target(*It, *g)].component) {
			auto e = boost::add_edge(source(*iterations.first, resultGraph), target(*iterations.first, resultGraph), *g);
		}
	}
	return resultGraph;
}

/**
* Tarjan's SCC algorithm
*/
void tarjanSCC(TarjanDirectedGraph* g) {

	auto v = new TarjanDirectedGraph::vertex_descriptor;
	*v = *boost::vertices(*g).first; //starting point of the algorithm, the first vertex in the graph
	tarjanStrongConnect(v, *g);
	tarjanI = 0;
	auto iterations = boost::vertices(*g);
	for(TarjanDirectedGraph g2 = *g; iterations.first < iterations.second; iterations.first++){
		if (g2[*iterations.first].number == -1) {
			TarjanDirectedGraph::vertex_descriptor v = *iterations.first;
			tarjanStrongConnect(&v, *g);
		}
	}

	/*
	 * memory deallocation
	 */
	
	delete v;

}

/**
* Tarjan's STRONGCONNECT procedure.
* implementation follows rigorously the pseudo-code
* input: a vertex, the graph, the stack, the auxiliary set, i, the component counter and the final result graph.
*/
void tarjanStrongConnect(TarjanDirectedGraph::vertex_descriptor* v, 
	TarjanDirectedGraph& g) {
	
	g[*v].lowpt = (tarjanI); //initialization as in the pseudo-code
	g[*v].lowvine = (tarjanI); //initialization as in the pseudo-code
	g[*v].number = (tarjanI); //initialization as in the pseudo-code
	tarjanI++;
	tarjanPointsStack.push(*v);
	tarjanPointsSet.insert(*v);

	auto iterations = adjacent_vertices(*v, g);
	auto It = iterations.first;

	for (; It < iterations.second; It++) {
		if (g[*It].number == -1) { // *It = w in the pseudo-code
			TarjanDirectedGraph::vertex_descriptor v2 = *It;
			tarjanStrongConnect(&v2, g);
			g[*v].lowpt = std::min(g[*v].lowpt, g[*It].lowpt);
			g[*v].lowvine = std::min(g[*v].lowvine, g[*It].lowvine);
		}
		// ELSE IF FOR ANCESTOR

		else if (g[*It].number<g[*v].number) {  // *It = w in the pseudo-code
			if (tarjanPointsSet.find(*It) != tarjanPointsSet.end()) {
				g[*v].lowvine = std::min(g[*v].lowvine, g[*It].number);
			}

		}
	}
	if ((g[*v].lowpt == g[*v].number) && (g[*v].lowvine == g[*v].number)) {
		tarjanComponentsCount++;

		std::cout << tarjanPointsStack.top() << "    " << g[tarjanPointsStack.top()].number << std::endl;
		std::cout << tarjanPointsStack.empty() << std::endl;

		while (!tarjanPointsStack.empty() && g[tarjanPointsStack.top()].number >= g[*v].number) {
			//SOMETHING WRONG
			g[tarjanPointsStack.top()].component = tarjanComponentsCount;
			tarjanPointsSet.erase(tarjanPointsStack.top());
			tarjanPointsStack.pop();
		}
		
	}
	
}


/**
* Nuutila's SCC algorithm
*/
void nuutilaSCC(NuutilaDirectedGraph* g) {

	std::cout << "vertices" << num_vertices(*g) << std::endl;
	nuutilaCounter = 0;
	for (auto iterations = boost::vertices(*g); iterations.first < iterations.second; iterations.first++) { // main loop of the algorithm
		NuutilaDirectedGraph::vertex_descriptor tempV = *iterations.first;
		nuutilaVisit(&tempV, *g);
	}	
	return;
}

/**
* Nuutila's VISIT2 procedure
* input: the vertex, the graph, the counter, the vertices Stack and the auxiliary vertices Set
* implementation follows rigorously the pseudo-code
*/
void nuutilaVisit(NuutilaDirectedGraph::vertex_descriptor* v,
	NuutilaDirectedGraph& g)
{

	g[*v].root = *v; // initialization as in the pseudo-code
	g[*v].inComponent = false; // initialization as in the pseudo-code
	g[*v].visited = true; // sets visited to true, to indicate the vertex has been visisted
	g[*v].visitIndex = nuutilaCounter; //sets the order of visist based on the counter
	nuutilaCounter++;


	auto iterations = adjacent_vertices(*v, g);
	auto It = iterations.first;

	//deque iterator not derefernciable
	while (It != iterations.second) { // while there are nodes adjacent to the initial node v
		if (g[*It].visited == false) { //visit them if they haven't been visited
			NuutilaDirectedGraph::vertex_descriptor tempV = *It;
			nuutilaVisit(&tempV, g);
		}
		if (!g[g[*It].root].inComponent) { // not inComponent[root[v]]
			g[g[*v].root].visitIndex = std::min(g[g[*v].root].visitIndex, g[g[*It].root].visitIndex);
		}
		It++;
	}
	assert(g[*v].root<5);
	if (g[*v].root == *v) { // root[v] == v
		if (!nuutilaVerticesStack.empty() && g[nuutilaVerticesStack.top()].visitIndex>g[*v].visitIndex) {
			assert(!nuutilaVerticesStack.empty());
			while (!nuutilaVerticesStack.empty() && g[nuutilaVerticesStack.top()].visitIndex < g[*v].visitIndex) {
				auto w = nuutilaVerticesStack.top();
				nuutilaVerticesSet.erase(w);
				assert(!nuutilaVerticesStack.empty());
				nuutilaVerticesStack.pop();
				g[w].inComponent = true;
			}
		}
		else {
			g[*v].inComponent = true;
		}
	}
	else if (nuutilaVerticesSet.find(g[*v].root) == nuutilaVerticesSet.end()) { // root[v] is not on the stack
		assert(g[*v].root<5);
		nuutilaVerticesStack.push(g[*v].root);
		nuutilaVerticesSet.insert(g[*v].root);
	}
	
	std::cout << nuutilaVerticesStack.top() << std::endl;
	return;
}

/**
* Imperative Version of Pearce's SCC algorithm. Based on PEA_FIND_SCC3
* implementation follows rigorously the pseudo-code
*/
void imperativePearceSCC(PearceDirectedGraph* g) {
	
	auto verticesStack = new std::stack<PearceDirectedGraph::vertex_descriptor>; //vS same as pseudo-code
	auto iteratorStack = new std::stack<int>; //iS same as pseudo-code
	auto c = new int; // c same as pseudo-code
	auto index = new int; //index same as pseudo-code

	*c = num_vertices(*g) - 1; // c initialization
	*index = 1; // index initialization

	for (auto iterations = boost::vertices(*g); iterations.first < iterations.second; iterations.first++) { // main loop of the algorithm
		PearceDirectedGraph::vertex_descriptor tempV = *iterations.first;
		pearceVisit(g, &tempV,verticesStack,iteratorStack,c,index);
	}

	//memory deallocation
	delete verticesStack;
	delete iteratorStack;
	delete c;
	delete index;
}

/**
* VISIT procedure
* input: the graph, the vertex, the vertices stack, the iterator stack, the int c, and the int index
* implementation follows rigorously the pseudo-code
*/
void pearceVisit(PearceDirectedGraph* g,
	PearceDirectedGraph::vertex_descriptor* v,
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, 
	std::stack<int>* iteratorStack,
	int* c,
	int* index) 
{
	pearceBeginVisiting(g, v, verticesStack, iteratorStack,c,index);
	while (!verticesStack->empty()) {
		pearceVisitLoop(g, verticesStack, iteratorStack, c, index);
	}
}

/** 
* BEGINVISIT procedure
* input: the graph, the vertex, the vertices stack, the iterator stack, the int c, and the int index
* implementation follows rigorously the pseudo-code
*/
void pearceBeginVisiting(PearceDirectedGraph* g, 
	PearceDirectedGraph::vertex_descriptor* v, 
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, 
	std::stack<int>* iteratorStack,
	int* c, 
	int* index) 
{
	verticesStack->push(*v);
	iteratorStack->push(0);

	PearceDirectedGraph tempGraph; // temporary auxiliary graph
	tempGraph = *g;
	tempGraph[*v].root = true;
	tempGraph[*v].rindex = *index;
	(*index)++;
	*g = tempGraph;

}

/**
* input: the graph, the vertices stack, the iterator stack, the int c, and the int index
* implementation follows rigorously the pseudo-code
*/
void pearceVisitLoop(PearceDirectedGraph* g,
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack,
	std::stack<int>* iteratorStack, 
	int* c, 
	int* index)
{
	auto v = new PearceDirectedGraph::vertex_descriptor;
	auto i = new int;
	*v = verticesStack->top();
	*i = iteratorStack->top();
	auto edgeLenght = new int;
	*edgeLenght = out_degree(*v, *g); // number of edges leaving the vertex
	while (*i <= *edgeLenght) { 
		if (*i > 0) {
			int j = *i - 1; // temporary variable
			pearceFinishEdge(g, v, &j);
		}
		if (*i < *edgeLenght && pearceBeginEdge(g, v, i, verticesStack, iteratorStack, c, index)) {
			return;
		}
		(*i)++;
	}

	pearceFinishVisiting(g, v, verticesStack, iteratorStack, c, index);
	
	// memory deallocation
	delete v;
	delete i;
	delete edgeLenght;
}

/**
* FINISHVISITING procedure
* input: the graph, the vertex, the vertices stack, the iterator stack, the int c, and the int index
* implementation follows rigorously the pseudo-code
*/
void pearceFinishVisiting(PearceDirectedGraph* g,
	PearceDirectedGraph::vertex_descriptor* v,
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, 
	std::stack<int>* iteratorStack, 
	int* c, 
	int* index) 
{
	verticesStack->pop();
	iteratorStack->pop();
	PearceDirectedGraph tempGraph; // temporary auxiliary graph
	tempGraph = *g;
	if (tempGraph[*v].root) { // root is true
		(*index)--;
		while (not(verticesStack->empty()) && tempGraph[*v].rindex <= tempGraph[verticesStack->top()].rindex) { //stack not empty and rindex[v] <= rindex[top(vS)]
			auto w = verticesStack->top(); // w in the pseudo-code
			verticesStack->pop();
			tempGraph[w].rindex = *c;
			(*index)--;
			*g = tempGraph;
		}
		tempGraph[*v].rindex = *c;
		(*c)--;
		*g = tempGraph;
	}
	else {
		verticesStack->push(*v);
	}
}

/**
* BEGINEDGE procedure
* input: the graph, the vertex, the index of the current outgoing edge i, the vertices stack, the iterator stack, the int c, and the int index
* implementation follows rigorously the pseudo-code
*/
bool pearceBeginEdge(PearceDirectedGraph* g, 
	PearceDirectedGraph::vertex_descriptor* v,
	int* i, 
	std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, 
	std::stack<int>* iteratorStack, 
	int* c, 
	int* index) 
{
	/*
	boost::graph_traits<PearceDirectedGraph>::adjacency_iterator ai;
	boost::graph_traits<PearceDirectedGraph>::adjacency_iterator ai_end;
	int j = 0;
	for (boost::tie(ai, ai_end) = adjacent_vertices(*v, *g); ai != ai_end, j != *i; ++ai, ++j){}
	*/
	auto iterations = adjacent_vertices(*v, *g);
	(int)*iterations.first;
	while (iterations.first < iterations.second && *iterations.first != *i) { //finds the vertex connected to the i-th edge of the node
		*iterations.first++; 
	}
	
	/*
	boost::graph_traits<PearceDirectedGraph>::vertex_iterator vi, vi_end, next;
	int j = 0;
	for (next = vi; vi != vi_end && j != *i; vi = next, j++) { //finds the vertex connected to the i-th edge of the node
	}
	*/
	

	PearceDirectedGraph tempGraph;
	tempGraph = *g;
	if (tempGraph[*iterations.first].rindex == 0) { // *iterations.first is w in the pseudo-code
		iteratorStack->pop();
		iteratorStack->push(*i + 1);
		auto tempV = *iterations.first;
		pearceBeginVisiting(g, &tempV, verticesStack, iteratorStack, c, index);
		*g = tempGraph;
		return true;
	}
	else {
		return false;
	}
}

/**
* FINISHEDGE procedure
* input: the graph, the vertex, the index of the current outgoing edge i
* implementation follows rigorously the pseudo-code
*/
void pearceFinishEdge(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, int* i) 
{
	PearceDirectedGraph tempGraph; // temporary auxiliary graph
	tempGraph = *g;

	auto iterations = adjacent_vertices(*v, *g);
	(int)*iterations.first;
	while (*iterations.first != *i && iterations.first < iterations.second) { //finds the vertex connected to the i-th edge of the node
		*iterations.first++;
	}
	if (tempGraph[*iterations.first].rindex<tempGraph[*v].rindex) { // *iterations.first is w in the pseudocode.
		tempGraph[*v].rindex = tempGraph[*iterations.first].rindex;
		tempGraph[*v].root = false;
		*g = tempGraph;
	}
}




