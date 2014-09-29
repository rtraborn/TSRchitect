#include "core/unused_api.h"
#include "extended/node_visitor_api.h"
#include "CagepNodeVisitor.h"

#define cagep_node_visitor_cast(GV) \
        gt_node_visitor_cast(cagep_node_visitor_class(), GV)


struct CagepNodeVisitor
{
  const GtNodeVisitor parent_instance;
  GtFeatureIndex *index;
  unsigned int *ngenes;
  unsigned int *nseqs;
};


/**
 * Cast a node visitor object as a CagepNodeVisitor
 *
 * @returns    a node visitor object cast as a CagepNodeVisitor
 */
const GtNodeVisitorClass* cagep_node_visitor_class();

/**
 * Destructor for the CagepNodeVisitor class
 *
 * @param[in] nv    the node visitor object
 */
static void cagep_node_visitor_free(GtNodeVisitor *nv);

/**
 * Validate and store any feature nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  fn       node representing a GFF3 feature entry
 * @param[out] error    error object to which error messages, if any, will be
 *                      written
 * @returns             0 for success, 1 for error
 */
static int cagep_node_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GT_UNUSED GtError *error);

/**
 * Store any region nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  rn       node representing a GFF3 sequence-region entry
 * @param[out] error    error object to which error messages, if any, will be
 *                      written
 * @returns             0 for success, 1 for error
 */
static int cagep_node_visitor_visit_region_node(GtNodeVisitor *nv,
                                                GtRegionNode *rn,
                                                GT_UNUSED GtError *error);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
const GtNodeVisitorClass* cagep_node_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (CagepNodeVisitor),
                                    cagep_node_visitor_free, NULL,
                                    cagep_node_visitor_visit_feature_node,
                                    cagep_node_visitor_visit_region_node, NULL,
                                    NULL);
  }
  return nvc;
}

static void cagep_node_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED CagepNodeVisitor *aiv = cagep_node_visitor_cast(nv);
}

GtNodeVisitor* cagep_node_visitor_new(GtFeatureIndex *index,
                                      unsigned int *ngenes, unsigned int *nseqs)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(cagep_node_visitor_class());
  CagepNodeVisitor *v = cagep_node_visitor_cast(nv);
  v->index = index;
  v->ngenes = ngenes;
  v->nseqs = nseqs;
  return nv;
}

static int cagep_node_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GT_UNUSED GtError *error)
{
  CagepNodeVisitor *v;
  gt_error_check(error);
  v = cagep_node_visitor_cast(nv);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *feat;
  for(feat = gt_feature_node_iterator_next(iter);
      feat != NULL;
      feat = gt_feature_node_iterator_next(iter))
  {
    if(gt_feature_node_has_type(feat, "gene"))
    {
      bool adderror = gt_feature_index_add_feature_node(v->index, feat, error);
      if(adderror)
        return 1;
      else
        *v->ngenes += 1;
    }
  }
  gt_feature_node_iterator_delete(iter);

  return 0;
}

static int cagep_node_visitor_visit_region_node(GtNodeVisitor *nv,
                                                GtRegionNode *rn,
                                                GT_UNUSED GtError *error)
{
  CagepNodeVisitor *v;
  gt_error_check(error);
  v = cagep_node_visitor_cast(nv);

  bool adderror = gt_feature_index_add_region_node(v->index, rn, error);
  if(adderror)
    return 1;
  else
    *v->nseqs += 1;

  return 0;
}
