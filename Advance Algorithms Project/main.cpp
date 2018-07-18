#include <iostream>
#include "boost\graph\adjacency_list.hpp"

boost::adjacency_list<>* createGraph(int numVertex);
std::unique_ptr<boost::adjacency_list<>> createGraph2(int numVertex);

int main() {
	std::cout << "Advanced Algorithms Project" << std::endl;
	
	boost::adjacency_list<> *g;
	g = createGraph(5);
	
	std::cin.get();
	return 0;
}



boost::adjacency_list<>* createGraph(int numVertex) {


	/*boost::adjacency_list<>* */ auto g = new boost::adjacency_list<>;
	/* boost::adjacency_list<>::vertex_descriptor* */ auto v = new boost::adjacency_list<>::vertex_descriptor[numVertex];
	for (int i = 0; i < numVertex; i++) {
		v[i] = add_vertex(*g);
	}
	for (int i = 0; i < numVertex; i++) {
		for (int j = i + 1; j < numVertex; j++) {
			auto e = boost::add_edge(i, j, *g);
		}
	}
	delete[] v;
	return g;

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






