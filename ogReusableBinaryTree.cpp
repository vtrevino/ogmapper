/* 
 * File:   ogReusableBinaryTree.cpp
 * Author: victortrevino
 * 
 * Created on July 12, 2022, 5:03 PM
 */

#include "ogReusableBinaryTree.hpp"

ogReusableBinaryTree::ogReusableBinaryTree() {
    nAllocated = 1000;
    nUsed = 0;
    blockSize = 100;
    nodes = NULL;
    keyDistance = 0;
    rootNode = -1;
    allocate();
}


ogReusableBinaryTree::ogReusableBinaryTree(uint32_t initialSize, uint32_t growingBlockSize) {
    nAllocated = initialSize;
    nUsed = 0;
    blockSize = growingBlockSize;
    nodes = NULL;
    keyDistance = 0;
    rootNode = -1;
    allocate();
}


ogReusableBinaryTree::~ogReusableBinaryTree() {
    //fprintf(stderr, "<Deallocating ogTree:"); fflush(stderr);
    //fprintf(stderr, "<~ogT:"); fflush(stderr);
    if (nodes != NULL) free(nodes);
    //fprintf(stderr, ":ogT>");
}


void ogReusableBinaryTree::allocate() {
    nodes = (ogrbtNode *) realloc(nodes, sizeof(ogrbtNODE) * nAllocated);
}



void ogReusableBinaryTree::allocateIfNeeded(uint32_t iNode) {
    if (iNode >= nAllocated) {
        nAllocated = iNode + blockSize;
        if (nAllocated / blockSize > 20) {
            blockSize = nAllocated / 10;
        }
        //fprintf(stderr, " %u ", nAllocated);
        allocate();
    }
}


void ogReusableBinaryTree::clear() {
    nUsed = 0;
    rootNode = -1;
}

void ogReusableBinaryTree::setMaximumDistanceForKey(uint32_t dist) {
    keyDistance = dist;
}

void ogReusableBinaryTree::setNodeTo(int32_t iNode, int64_t theKey, int64_t theContent) {
    allocateIfNeeded(iNode);
    ogrbtNode *p = nodes+iNode;
    p->key = theKey;
    p->content = theContent;
    p->left = -1;
    p->right = -1;
}


void ogReusableBinaryTree::setNodeContent(int32_t iNode, int64_t theContent) {
    allocateIfNeeded(iNode);
    nodes[iNode].content = theContent;
}

int64_t ogReusableBinaryTree::incNodeContent(int32_t iNode) {
    allocateIfNeeded(iNode);
    return ++nodes[iNode].content;
}

int64_t ogReusableBinaryTree::getNodeContent(int32_t iNode) {
    if (iNode < nAllocated) {
        return nodes[iNode].content;
    } else {
        return 0;
    }
}

int64_t ogReusableBinaryTree::getNodeKey(int32_t iNode) {
    if (iNode < nAllocated) {
        return nodes[iNode].key;
    } else {
        return 0;
    }
}

int32_t ogReusableBinaryTree::searchForKey(int64_t theKey) {
    if (nUsed == 0) return -1; 
    int32_t i = rootNode;
    int32_t last_i = 0;
    int64_t d;
    ogrbtNode *p;
    while (i >= 0) {
        last_i = i;
        p = nodes+i;
        d = p->key - theKey;
        if (d < 0) d = -d;
        if (d <= keyDistance) {
            return i;
        }
        i = (p->key > theKey ? p->left : p->right);
    }
    return -last_i-1;
}

int32_t ogReusableBinaryTree::insertNode(int64_t theKey, int64_t defaultContent) {
    int32_t i = searchForKey(theKey);
    if (i >= 0) {
        return i;
    }
    allocateIfNeeded(nUsed+1);
    ogrbtNode *p = nodes + nUsed;
    p->key = theKey;
    p->content = defaultContent;
    p->left = -1;
    p->right = -1;
    if (nUsed > 0) {
        i = (-i)-1;
        if (nodes[i].key > theKey) {
            if (nodes[i].left != -1) {
                fprintf(stderr, "***** PROBLEM in node[i].left i=%d\n",i);
            }
            nodes[i].left = nUsed;
        } else {
            if (nodes[i].right != -1) {
                fprintf(stderr, "***** PROBLEM in node[i].right i=%d\n",i);
            }
            nodes[i].right = nUsed;
        }
        i = -nUsed;
    } else {
        rootNode = 0;
        i = 0;
    }
    nUsed++;
    return i;
}


void ogReusableBinaryTree::buildTreeFromUnorderedNodesSet(uint32_t nNodes) {
    allocateIfNeeded(nNodes);
    nUsed = 0;
    rootNode = -1;
    uint32_t i;
    for (i = 0; i < nNodes; i++) {
        insertNode(nodes[i].key, nodes[i].content);
    }
}

void ogReusableBinaryTree::buildTreeFromOrderedNodesSet(uint32_t nNodes) {
    nUsed = nNodes;
    if (nUsed > 2) {
        rootNode = nUsed >> 1;
        setNodes(rootNode, 0, rootNode-1, rootNode+1, nUsed-1);
    } else if (nUsed == 2) {
        rootNode = 1;
        nodes[1].left = 0;
        nodes[1].right = -1;
    } else if (nUsed == 1) {
        rootNode = 0;
        nodes[0].left = -1;
        nodes[0].right = -1;
    } else {
        rootNode = -1;
    }
}

void ogReusableBinaryTree::setNodes(int32_t parent, int32_t leftFrom, int32_t leftTo, int32_t rightFrom, int32_t rightTo) {
    allocateIfNeeded(parent);
    //fprintf(stderr, "*****\n");
    //fprintf(stderr, "Parent %d, Left From %d, Left To %d, Right From %d, Right To %d\n", parent, leftFrom, leftTo, rightFrom, rightTo);    
    int32_t left = (leftFrom + leftTo) >> 1;
    int32_t right = (rightFrom + rightTo) >> 1;
    nodes[parent].left = left;
    nodes[parent].right = right;
    if (left > leftFrom && left < leftTo) {
        setNodes(left, leftFrom, left - 1, left + 1, leftTo);
    } else if (leftTo == leftFrom) {
        nodes[left].left = -1;
        nodes[left].right = -1;
    } else {
        // must: left == leftFrom
        nodes[left].left = -1;
        nodes[left].right = leftTo;
    }
    if (right > rightFrom && right < rightTo) {
        setNodes(right, rightFrom, right - 1, right + 1, rightTo);
    } else if (rightTo == rightFrom) {
        nodes[right].left = -1;
        nodes[right].right = -1;
    } else {
        // must: right == rightFrom
        nodes[right].left = -1;
        nodes[right].right = rightTo;
    }
}

void ogReusableBinaryTree::printNodesInContentRange(uint64_t minValue, uint64_t maxValue) {
    printNodesInContentRange(rootNode, minValue, maxValue);
}

void ogReusableBinaryTree::printNodesInContentRange(int32_t iNode, uint64_t minValue, uint64_t maxValue) {
    if (nodes[iNode].left >= 0) {
        printNodesInContentRange(nodes[iNode].left, minValue, maxValue);
    }
    if (nodes[iNode].content >= minValue && nodes[iNode].content <= maxValue) {
        fprintf(stderr, "Node %d, Key %lld, Content %lld\n", iNode, nodes[iNode].key, nodes[iNode].content);
    }
    if (nodes[iNode].right >= 0) {
        printNodesInContentRange(nodes[iNode].right, minValue, maxValue);
    }
}