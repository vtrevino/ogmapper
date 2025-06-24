/* 
 * File:   ogReusableBinaryTree.hpp
 * Author: victortrevino
 *
 * Created on July 12, 2022, 5:03 PM
 */

#ifndef OGREUSABLEBINARYTREE_HPP
#define OGREUSABLEBINARYTREE_HPP

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct ogrbtNODE {
    int64_t     key;
    int64_t     content;
    int32_t     left, right; // index positions relative to root, to avoid redoing in realloc
} ogrbtNode;

class ogReusableBinaryTree {
    uint32_t    blockSize;
    uint32_t    keyDistance;
    int32_t     rootNode;
    ogrbtNode   *nodes;
    
    
public:
    uint32_t    nAllocated;
    uint32_t    nUsed;
                            ogReusableBinaryTree();
                            ogReusableBinaryTree(uint32_t initialSize, uint32_t growingBlockSize);
    virtual                ~ogReusableBinaryTree();
    void                    setMaximumDistanceForKey(uint32_t dist);
    void                    clear(); // remove all nodes
    void                    allocate();
    void                    allocateIfNeeded(uint32_t iNode);
    void                    setNodeTo(int32_t iNode, int64_t theKey, int64_t theContent);
    void                    setNodeContent(int32_t iNode, int64_t theContent);
    int64_t                 incNodeContent(int32_t iNode);
    int64_t                 getNodeContent(int32_t iNode);
    int64_t                 getNodeKey(int32_t iNode);
    int32_t                 searchForKey(int64_t theKey); // get Node index or -1 if not found
    int32_t                 insertNode(int64_t theKey, int64_t defaultContent);
    void                    buildTreeFromOrderedNodesSet(uint32_t nNodes);
    void                    setNodes(int32_t parent, int32_t leftFrom, int32_t leftTo, int32_t rightFrom, int32_t rightTo);
    void                    buildTreeFromUnorderedNodesSet(uint32_t nNodes);
    void                    printNodesInContentRange(uint64_t minValue, uint64_t maxValue);
    void                    printNodesInContentRange(int32_t iNode, uint64_t minValue, uint64_t maxValue);
    
private:

};

#endif /* OGREUSABLEBINARYTREE_HPP */

