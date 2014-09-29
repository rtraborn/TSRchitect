#ifndef CAGEP_NODE_VISITOR
#define CAGEP_NODE_VISITOR

#include "genometools.h"

typedef struct CagepNodeVisitor CagepNodeVisitor;

/**
 * Allocate memory for a genome node visitor used to load gene annotations into
 * memory.
 *
 * @param[out] index     the feature index into which data will be loaded
 * @param[out] ngenes    counter for the number of gene annotations loaded into
 *                       memory
 * @param[out] nseqs     count for the number of sequences annotated by the
 *                       genes loaded into memory
 * @returns              a node visitor object
 */
GtNodeVisitor* cagep_node_visitor_new(GtFeatureIndex *index,
                                     unsigned int *ngenes, unsigned int *nseqs);

#endif
