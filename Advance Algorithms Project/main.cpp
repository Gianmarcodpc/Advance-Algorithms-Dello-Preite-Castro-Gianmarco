#include <iostream>
#include <boost\graph\adjacency_list.hpp>
#include <time.h>
#include <boost\random\mersenne_twister.hpp>
#include <boost\random\uniform_int_distribution.hpp>
#include <boost\graph\graphviz.hpp>
#include <string>

typedef boost::adjacency_list	<boost::vecS, boost::vecS,boost::directedS,boost::no_property,boost::no_property> DirectedGraph;


DirectedGraph* createRandomDirectedGraph(int numVertex);
std::unique_ptr<boost::adjacency_list<>> createGraph2(int numVertex);
void saveDirectedGraphToFile(DirectedGraph* g,int numVertex, int index);
void generateAndSaveRandomGraphs(int numGraphs);


int main() {
	std::cout << "Advanced Algorithms Project" << std::endl;
	
	generateAndSaveRandomGraphs(10);
	std::cin.get();
	return 0;
}


DirectedGraph* createRandomDirectedGraph(int numVertex) {

	boost::random::mt19937 gen(time(NULL));
	boost::random::uniform_int_distribution<> dist(1, 10000);
	int treshold = dist(gen);

	/*DirectedGraph* */ auto g = new DirectedGraph;
	/* boost::adjacency_list<>::vertex_descriptor* */ auto v = new boost::adjacency_list<>::vertex_descriptor[numVertex];
	for (int i = 0; i < numVertex; i++) {
		v[i] = add_vertex(*g);
	}
	for (int i = 0, randValue = 0; i < numVertex; i++) {
		for (int j = 0; j < numVertex; j++) {
			randValue = dist(gen);
			if (randValue >= treshold) {
				auto e = boost::add_edge(i, j, *g);
			}
		}
	}
	delete[] v;
	return g;

}

void saveDirectedGraphToFile(DirectedGraph* g, int numVertex, int index) {
	

	std::string fileName = "DirectedGraph" + std::to_string(numVertex) + "Vertexes" + std::to_string(index);
	std::string path = "./DirectedGraphs/" + fileName + ".txt";
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


std::unique_ptr<boost::adjacency_list<>> createGraph2(int numVertex) {

	std::unique_ptr<boost::adjacency_list<>> g = std::make_unique<boost::adjacency_list<>>();
	/* boost::adjacency_list<>::vertex_descriptor* */ auto v = new boost::adjacency_list<>::vertex_descriptor[numVertex];
	for (int i = 0; i < numVertex; i++) {
		v[i] = add_vertex(*g);
	}
	delete[] v;

	return g;

}






