#ifndef TREE_H_
#define TREE_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include "Helper.h"
#include "Node.h"


#ifndef NULL
#define NULL 0
#endif


/*
 * Source (Def. & Decl.)
 */

template <class T, class C >
class Tree {
public:
	Tree() : root(NULL) {};
	Tree(Node<T, C> *root) : root(root) {};
	virtual ~Tree() {
		if (this->root != NULL) delete this->root;
	};

	inline void setRoot(Node<T, C> *root) {this->root = root;};
	inline Node<T, C> *getRoot() const {return this->root;};


	virtual T classify(const std::vector<T > &sample) const {
		return this->root->classify(sample);
	}

	virtual std::ostream& print(std::ostream& os) const {
		this->root->print(os);
		return os;
	}

	virtual bool isThereVarID(uli_t varID) const {
    return this->root->isThereVarID(varID);
	}

	virtual void printXml(std::ostream& os) const {
    os
    << "<node id=\"0\" "
    << "size=\"" << this->root->size() << "\">" << std::endl;
		this->root->printXml(os);
		os << "</node>" << std::endl;
	}

	virtual void summary() const {
		std::ostream& os = std::cout;
//		os	<< "Tree:" << std::endl;
//		this->root->print(os);

		os	<< "Number of leafs: "
			<< this->root->getNumOfLeafs() << std::endl;
/*
		os	<< "Max depth: "
			<< this->root->getMaxDepth() << std::endl;
*/
		os	<< "Number of nodes: "
			<< this->root->getNumOfNodes() << std::endl;
        }

	Node<T, C> *root;
};

template <class T, class C>
std::ostream &operator<<(std::ostream &os, Tree<T, C > &tree) {
	return tree.print(os);
}

#endif /*TREE_H_*/
