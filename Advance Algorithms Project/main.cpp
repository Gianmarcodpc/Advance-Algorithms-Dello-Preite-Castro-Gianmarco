#include <iostream>
#include <boost\graph\adjacency_list.hpp>
#include <boost\graph\adjacency_iterator.hpp>
#include <time.h>
#include <boost\random\mersenne_twister.hpp>
#include <boost\random\uniform_int_distribution.hpp>
#include <boost\graph\graphviz.hpp>
#include <string>
#include <stack>
#include <set>



struct TarjanVertexData {
	int lowpt;
	int lowvine;
	int number = -1;
	int component;
};

struct NuutilaVertexData {
	boost::adjacency_list<>::vertex_descriptor root;
	bool inComponent;
	bool visited = false;
	int visitIndex;
};

struct PearceVertexData {
	int rindex;
	bool root = false;
};



typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::no_property> DirectedGraph;
typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, TarjanVertexData, boost::no_property> TarjanDirectedGraph;
typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, NuutilaVertexData, boost::no_property> NuutilaDirectedGraph;
typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, PearceVertexData, boost::no_property> PearceDirectedGraph;


DirectedGraph* createRandomDirectedGraph(int numVertex);
TarjanDirectedGraph* createRandomTarjanDirectedGraph(int numVertex);
NuutilaDirectedGraph* createRandomNuutilaDirectedGraph(int numVertex);
PearceDirectedGraph* createRandomPearceDirectedGraph(int numVertex);
void saveDirectedGraphToFile(DirectedGraph* g,int numVertex, int index);
void saveTarjanDirectedGraphToFile(TarjanDirectedGraph* g, int numVertex, int index);
void generateAndSaveRandomGraphs(int numGraphs);
TarjanDirectedGraph* tarjanSCC(TarjanDirectedGraph* g);
void tarjanStrongConnect(boost::adjacency_list<>::vertex_descriptor* v, TarjanDirectedGraph* g, std::stack<TarjanDirectedGraph::vertex_descriptor>* points, std::set<TarjanDirectedGraph::vertex_descriptor>* pointsSet, int* i,int* componentCounter, TarjanDirectedGraph* resultGraph);
TarjanDirectedGraph* displayTarjanSCC(TarjanDirectedGraph* g);


void nuutilaSCC(NuutilaDirectedGraph* g);
void nuutilaVisit(NuutilaDirectedGraph::vertex_descriptor* v, NuutilaDirectedGraph* g, int* counter, std::stack<NuutilaDirectedGraph::vertex_descriptor>* verticesStack, std::set<NuutilaDirectedGraph::vertex_descriptor>* verticesSet);


void imperativePearceSCC(PearceDirectedGraph* g);
void pearceVisit(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index);
void pearceBeginVisiting(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index);
void pearceVisitLoop(PearceDirectedGraph* g, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index);
void pearceFinishVisiting(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index);
bool pearceBeginEdge(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, int* i, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index);
void pearceFinishEdge(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, int* i);





int main() {
	std::cout << "Advanced Algorithms Project" << std::endl;



	PearceDirectedGraph* g = createRandomPearceDirectedGraph(5);
	imperativePearceSCC(g);
	//NuutilaDirectedGraph *g = createRandomNuutilaDirectedGraph(5);
	//nuutilaSCC(g);
	//TarjanDirectedGraph* g = createRandomTarjanDirectedGraph(5);
	//TarjanDirectedGraph* sccGraph = displayTarjanSCC(g);
	//delete sccGraph;
	//generateAndSaveRandomGraphs(10);
	std::cin.get();
	return 0;
}


DirectedGraph* createRandomDirectedGraph(int numVertex) {

	boost::random::mt19937 gen(time(NULL));
	boost::random::uniform_int_distribution<> dist(1, 10000);
	int treshold = dist(gen);
	
	/*DirectedGraph* */ auto g = new DirectedGraph;
	/* boost::adjacency_list<>::vertex_descriptor* */ auto vertices = new boost::adjacency_list<>::vertex_descriptor[numVertex];
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

TarjanDirectedGraph* createRandomTarjanDirectedGraph(int numVertex) {

	boost::random::mt19937 gen(1);
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

NuutilaDirectedGraph* createRandomNuutilaDirectedGraph(int numVertex) {

	boost::random::mt19937 gen(1);
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

void saveDirectedGraphToFile(DirectedGraph* g, int numVertex, int index) {


	std::string fileName = "DirectedGraph" + std::to_string(numVertex) + "Vertexes" + std::to_string(index);
	std::string path = "./DirectedGraphs/" + fileName + ".txt";
	std::ofstream writer(path);
	boost::write_graphviz(writer, *g);
	return;
}

void saveTarjanDirectedGraphToFile(TarjanDirectedGraph* g, int numVertex, int index) {
	

	std::string fileName = "TarjanDirectedGraph" + std::to_string(numVertex) + "Vertexes" + std::to_string(index);
	std::string path = "./TarjanDirectedGraphs/" + fileName + ".txt";
	std::ofstream writer(path);
	boost::write_graphviz(writer, *g);
	return;
}

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

TarjanDirectedGraph* displayTarjanSCC(TarjanDirectedGraph* g) {
	TarjanDirectedGraph* resultGraph = tarjanSCC(g);
	TarjanDirectedGraph tempGraph;
	tempGraph = *resultGraph;
	for (auto iterations = boost::edges(*g); iterations.first != iterations.second; iterations.first++) {
		if (tempGraph[source(*iterations.first, *g)].component == tempGraph[target(*iterations.first, *g)].component) {
			auto e = boost::add_edge(source(*iterations.first, *g), target(*iterations.first, *g), *g);
		}
	}
	*resultGraph = tempGraph;
	return resultGraph;
}

TarjanDirectedGraph* tarjanSCC(TarjanDirectedGraph* g) {
	auto i = new int;
	*i = 0;
	auto componentCounter = new int;
	*componentCounter = 0;
	auto v = new TarjanDirectedGraph::vertex_descriptor;
	*v = *boost::vertices(*g).first;
	auto points = new std::stack<TarjanDirectedGraph::vertex_descriptor>;
	auto pointsSet = new std::set<TarjanDirectedGraph::vertex_descriptor>;
	auto resultGraph = new TarjanDirectedGraph;
	tarjanStrongConnect(v, g, points, pointsSet, i, componentCounter, resultGraph);
	*i = 0;
	auto iterations = boost::vertices(*g);
	for(TarjanDirectedGraph g2 = *g; iterations.first < iterations.second; iterations.first++){
		if (g2[*iterations.first].number == -1) {
			*g = g2;
			TarjanDirectedGraph::vertex_descriptor v = *iterations.first;
			tarjanStrongConnect(&v, g, points, pointsSet, i, componentCounter, resultGraph);
		}
	}

	delete i;
	delete componentCounter;
	delete v;
	delete points;
	delete pointsSet;
	return resultGraph;

}


void tarjanStrongConnect(TarjanDirectedGraph::vertex_descriptor* v, TarjanDirectedGraph* g, std::stack<TarjanDirectedGraph::vertex_descriptor>* points, std::set<TarjanDirectedGraph::vertex_descriptor>* pointsSet, int* i,int* componentCounter, TarjanDirectedGraph* resultGraph) {
	TarjanDirectedGraph g2 = *g;
	g2[*v].lowpt = (*i);
	g2[*v].lowvine = (*i);
	g2[*v].number = (*i);
	(*i)++;
	*g = g2;
	points->push(*v);
	pointsSet->insert(*v);
	TarjanDirectedGraph::vertex_descriptor top;

	
	for (auto iterations = adjacent_vertices(*v, *g); iterations.first < iterations.second; iterations.first++) {
		if (g2[*iterations.first].number == -1) {
			TarjanDirectedGraph::vertex_descriptor v2 = *iterations.first;
			tarjanStrongConnect(&v2, g, points, pointsSet, i, componentCounter, resultGraph);
			g2[*v].lowpt = std::min(g2[*v].lowpt, g2[*iterations.first].lowpt);
			g2[*v].lowvine = std::min(g2[*v].lowvine, g2[*iterations.first].lowvine);
			*g = g2;
		}
		// ELSE IF FOR ANCESTOR

		else if (g2[*iterations.first].number<g2[*v].number) {
			if (pointsSet->find(*iterations.first) == pointsSet->end()) {
				g2[*v].lowvine = std::min(g2[*v].lowvine, g2[*iterations.first].number);
				*g = g2;
			}

		}
	}
	if ((g2[*v].lowpt == g2[*v].number)&&(g2[*v].lowvine== g2[*v].number)) {
		(*componentCounter)++;
		top = points->top();
		while (g2[top].number>=g2[*v].number) {
			//SOMETHING WRONG
			resultGraph->added_vertex(top);
			TarjanDirectedGraph g3 = *resultGraph;
			g3[top].component = *componentCounter;
			*resultGraph = g3;
			pointsSet->erase(top);
			points->pop();
			top = points->top();
		}
	}
	
}


void nuutilaSCC(NuutilaDirectedGraph* g) {
	auto counter = new int;
	*counter = 0;
	auto verticesStack = new std::stack<NuutilaDirectedGraph::vertex_descriptor>;
	auto verticesSet = new std::set<NuutilaDirectedGraph::vertex_descriptor>;
	NuutilaDirectedGraph fakeGraph;
	auto fakeVertex = add_vertex(fakeGraph);
	fakeGraph[fakeVertex].visitIndex = 0;
	verticesStack->push(fakeVertex);
	verticesSet->insert(fakeVertex);
	
	for (auto iterations = boost::vertices(*g); iterations.first < iterations.second; iterations.first++) {
		NuutilaDirectedGraph::vertex_descriptor tempV = *iterations.first;
		nuutilaVisit(&tempV, g, counter, verticesStack, verticesSet);
	}
	//add all the deletes
	delete counter;
	delete verticesStack;
	delete verticesSet;


}

void nuutilaVisit(NuutilaDirectedGraph::vertex_descriptor* v, NuutilaDirectedGraph* g, int* counter, std::stack<NuutilaDirectedGraph::vertex_descriptor>* verticesStack, std::set<NuutilaDirectedGraph::vertex_descriptor>* verticesSet) {
	NuutilaDirectedGraph tempGraph = *g;
	tempGraph[*v].root = *v;
	tempGraph[*v].inComponent = false; 
	tempGraph[*v].visited = true;
	tempGraph[*v].visitIndex = *counter;
	(*counter)++;
	verticesStack->push(*v);
	verticesSet->insert(*v);
	*g = tempGraph;
	auto iterations = adjacent_vertices(*v, tempGraph);
	//THERE IS AN ERROR HERE. DEREFERENCING ITERATOR
	while(iterations.first != iterations.second){	
		if (tempGraph[*iterations.first].visited == false) {
			NuutilaDirectedGraph::vertex_descriptor tempV = *iterations.first;
			nuutilaVisit(&tempV, &tempGraph, counter, verticesStack, verticesSet); // *g instead of tempGraph?
		}
		//auto tempV2 = tempGraph[*iterations.first].root;
		if (!tempGraph[tempGraph[*iterations.first].root].inComponent) {
			tempGraph[*v].root = std::min(tempGraph[*v].visitIndex, tempGraph[*iterations.first].visitIndex);
		}
		if (tempGraph[*v].root == *v) {
			if (tempGraph[verticesStack->top()].visitIndex>tempGraph[*v].visitIndex) {
				while (tempGraph[verticesStack->top()].visitIndex >= tempGraph[*v].visitIndex) {
					auto w = verticesStack->top();
					verticesSet->erase(w);
					verticesStack->pop();
					tempGraph[w].inComponent = true;
				}
			}
			else {
				tempGraph[*v].inComponent = true;
			}
		}
		else if(verticesSet->find(tempGraph[*v].root) == verticesSet->end()){
			verticesStack->push(tempGraph[*v].root);
			verticesSet->insert(tempGraph[*v].root);
		}
		iterations.first++;
	}
	*g = tempGraph;
}




void imperativePearceSCC(PearceDirectedGraph* g) {
	
	auto verticesStack = new std::stack<PearceDirectedGraph::vertex_descriptor>;
	auto iteratorStack = new std::stack<int>;
	auto c = new int;
	auto index = new int;

	*c = num_vertices(*g) - 1;
	*index = 1;

	for (auto iterations = boost::vertices(*g); iterations.first < iterations.second; iterations.first++) {
		PearceDirectedGraph::vertex_descriptor tempV = *iterations.first;
		pearceVisit(g, &tempV,verticesStack,iteratorStack,c,index);
	}

	delete verticesStack;
	delete iteratorStack;
	delete c;
	delete index;
}


void pearceVisit(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index) {
	pearceBeginVisiting(g, v, verticesStack, iteratorStack,c,index);
	while (!verticesStack->empty()) {
		pearceVisitLoop(g, verticesStack, iteratorStack, c, index);
	}
}

void pearceBeginVisiting(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack,int* c, int* index) {
	verticesStack->push(*v);
	iteratorStack->push(0);
	PearceDirectedGraph tempGraph;
	tempGraph = *g;
	tempGraph[*v].root = true;
	tempGraph[*v].rindex = *index;
	(*index)++;
	*g = tempGraph;

}

void pearceVisitLoop(PearceDirectedGraph* g, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index) {
	auto v = new PearceDirectedGraph::vertex_descriptor;
	auto i = new int;
	*v = verticesStack->top();
	*i = iteratorStack->top();
	auto edgeLenght = new int;
	*edgeLenght = out_degree(*v, *g);
	while (*i <= *edgeLenght) {
		if (*i > 0) {
			int j = *i - 1;
			pearceFinishEdge(g, v, &j);
		}
		if (*i < *edgeLenght && pearceBeginEdge(g, v, i, verticesStack, iteratorStack, c, index)) {
			return;
		}
		(*i)++;
	}

	pearceFinishVisiting(g, v, verticesStack, iteratorStack, c, index);
	delete v;
	delete i;
	delete edgeLenght;
}

void pearceFinishVisiting(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index) {
	verticesStack->pop();
	iteratorStack->pop();
	PearceDirectedGraph tempGraph;
	tempGraph = *g;
	if (tempGraph[*v].root) {
		(*index)--;
		while (not(verticesStack->empty()) && tempGraph[*v].rindex <= tempGraph[verticesStack->top()].rindex) {
			auto w = verticesStack->top();
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

bool pearceBeginEdge(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, int* i, std::stack<PearceDirectedGraph::vertex_descriptor>* verticesStack, std::stack<int>* iteratorStack, int* c, int* index) {
	/*
	auto iterations = adjacent_vertices(*v, *g);
	(int)*iterations.first;
	while (*iterations.first != *i && iterations.first < iterations.second) {
		*iterations.first++; 
	}
	*/
	
	boost::graph_traits<PearceDirectedGraph>::vertex_iterator vi, vi_end, next;
	int j = 0;
	for (next = vi; vi != vi_end && j != *i; vi = next, j++) {
	}
	

	PearceDirectedGraph tempGraph;
	tempGraph = *g;
	if (tempGraph[*vi].rindex == 0) {
		iteratorStack->pop();
		iteratorStack->push(*i + 1);
		auto tempV = *vi;
		pearceBeginVisiting(g, &tempV, verticesStack, iteratorStack, c, index);
		*g = tempGraph;
		return true;
	}
	else {
		return false;
	}
}

void pearceFinishEdge(PearceDirectedGraph* g, PearceDirectedGraph::vertex_descriptor* v, int* i) {
	PearceDirectedGraph tempGraph;
	tempGraph = *g;
	auto iterations = adjacent_vertices(*v, *g);
	for (; iterations.first < iterations.second && *iterations.first != (size_t)*i; iterations.first++) {

	}
	if (tempGraph[*iterations.first].rindex<tempGraph[*v].rindex) {
		tempGraph[*v].rindex = tempGraph[*iterations.first].rindex;
		tempGraph[*v].root = false;
		*g = tempGraph;
	}
}




