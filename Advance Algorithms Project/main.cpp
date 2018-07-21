#include <iostream>
#include <boost\graph\adjacency_list.hpp>
#include <time.h>
#include <boost\random\mersenne_twister.hpp>
#include <boost\random\uniform_int_distribution.hpp>
#include <boost\graph\graphviz.hpp>
#include <string>

typedef boost::adjacency_list	<boost::vecS, boost::vecS,boost::directedS,boost::no_property,boost::no_property> DirectedGraph;


DirectedGraph* createRandomGraph(int numVertex);
std::unique_ptr<boost::adjacency_list<>> createGraph2(int numVertex);
void saveDirectedGraphToFile(DirectedGraph* g,int numVertex);


int main() {
	std::cout << "Advanced Algorithms Project" << std::endl;
	DirectedGraph* g;
	g = createRandomGraph(5);
	saveDirectedGraphToFile(g, 5);

	delete g;
	std::cin.get();
	return 0;
}



DirectedGraph* createRandomGraph(int numVertex) {

	boost::random::mt19937 gen(time(NULL));
	boost::random::uniform_int_distribution<> dist(1, 10000);
	int treshold = dist(gen);
	std::cout << "Treshold " << treshold << std::endl;

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
				std::cout << "Added Edge from " << i << " to " << j << std::endl;
			}
			std::cout << "randValue " << j+i*numVertex << " " << randValue << std::endl;
		}
	}
	delete[] v;
	return g;

}

void saveDirectedGraphToFile(DirectedGraph* g, int numVertex) {
	
	std::ofstream writer("./DirectedGraphs/MyFirstDirectedGraph.txt");
	boost::write_graphviz(writer, *g);
	return;
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






