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


typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::no_property> DirectedGraph;
typedef boost::adjacency_list	<boost::vecS, boost::vecS, boost::directedS, TarjanVertexData, boost::no_property> TarjanDirectedGraph;


DirectedGraph* createRandomDirectedGraph(int numVertex);
TarjanDirectedGraph* createRandomTarjanDirectedGraph(int numVertex);
void saveDirectedGraphToFile(DirectedGraph* g,int numVertex, int index);
void saveTarjanDirectedGraphToFile(TarjanDirectedGraph* g, int numVertex, int index);
void generateAndSaveRandomGraphs(int numGraphs);
TarjanDirectedGraph* tarjanSCC(TarjanDirectedGraph* g);
void tarjanStrongConnect(boost::adjacency_list<>::vertex_descriptor* v, TarjanDirectedGraph* g, std::stack<TarjanDirectedGraph::vertex_descriptor>* points, std::set<TarjanDirectedGraph::vertex_descriptor>* pointsSet, int* i,int* componentCounter, TarjanDirectedGraph* resultGraph);
TarjanDirectedGraph* displayTarjanSCC(TarjanDirectedGraph* g);

int main() {
	std::cout << "Advanced Algorithms Project" << std::endl;


	TarjanDirectedGraph* g = createRandomTarjanDirectedGraph(5);
	TarjanDirectedGraph* sccGraph = displayTarjanSCC(g);
	delete sccGraph;
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

	boost::random::mt19937 gen(time(NULL));
	boost::random::uniform_int_distribution<> dist(1, 10000);
	int treshold = dist(gen);

	/*DirectedGraph* */ auto g = new TarjanDirectedGraph;
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
	auto i = new int();
	*i = 0;
	auto componentCounter = new int();
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
	auto iterations = adjacent_vertices(*v, *g);
	for (; iterations.first < iterations.second; iterations.first++) {
		if (g2[*iterations.first].number == -1) {
			tarjanStrongConnect(v, g, points, pointsSet, i, componentCounter, resultGraph);
			g2[*v].lowpt = std::min(g2[*v].lowpt, g2[*iterations.first].lowpt);
			g2[*v].lowvine = std::min(g2[*v].lowvine, g2[*iterations.first].lowvine);
			*g = g2;
		}
		else if (g2[*iterations.first].number<g2[*v].number) {
			if (pointsSet->find(*iterations.first) == pointsSet->end()) {
				g2[*v].lowvine = std::min(g2[*v].lowvine, g2[*iterations.first].number);
				*g = g2;
			}

		}
	}
	if ((g2[*v].lowpt == g2[*v].number)&&(g2[*v].lowvine== g2[*v].number)) {
		*componentCounter++;
		top = points->top();
		while (g2[top].number>= g2[*v].number) {
			resultGraph->added_vertex(points->top());
			TarjanDirectedGraph g3 = *resultGraph;
			g3[points->top()].component = *componentCounter;
			pointsSet->erase(points->top());
			points->pop();
		}
	}
	
}

