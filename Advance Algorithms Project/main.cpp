#include <iostream>
#include "boost\graph\adjacency_list.hpp"



int main() {
	std::cout << "Advanced Algorithms Project" << std::endl;
	std::cin.get();
	boost::adjacency_list<> g;
	auto v1 = add_vertex(g);
	auto v2 = add_vertex(g);
	boost::adjacency_list<>::vertex_descriptor v3 = add_vertex(g);
	
	
	
	return 0;
}


boost::adjacency_list<>* createGraph(int numVertex) {

	/*boost::adjacency_list<>* */ auto g = new boost::adjacency_list<>;
	/* boost::adjacency_list<>::vertex_descriptor* */ auto v = new boost::adjacency_list<>::vertex_descriptor[numVertex];
	for (int i = 0; i < numVertex; i++) {
		v[i] = add_vertex(*g);
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






