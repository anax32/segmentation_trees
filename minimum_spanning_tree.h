/**
 * MSTrees
 */
#define OP_QMAP  1

class MSTree
{
public:
  typedef struct edge_s
  {
    unsigned int v1, v2;
    float w;
  } edge_t;

  static bool edgeComp (edge_t e1, edge_t e2)  {return (e1.w<e2.w);}

protected:
  const unsigned int V; ///< number of vertices in the graph
  const unsigned int E; ///< number of edges in the full description
  const unsigned int *e; ///< edge list as a pair of vertex indicies
  const float *w; ///< edge weights as an array of floats
  const unsigned int TE; ///< total edges in the mst
  const unsigned int TC; ///< total components in the hierarchy
  Tree **mst; ///< tree describing the component map
  const bool verbose; ///< verbose output
  edge_t *edgelist; ///< sortable edgelist

  void allocTree ()
  {
    if (*mst == NULL)
    {
      *mst = new Tree ();
    }

    if ((*mst)->length () == 0)
    {
      (*mst)->setLength (TC);
    }
  }

  void initialComponents ()
  {
    unsigned int  i;

    // each component is intially its own parent (a forest)
    for (i=0; i<TC; i++)
    {
      (*mst)->get (i)->setCID (i);
    }
  }

  void allocEdgeList ()
  {
    unsigned int  i;

    edgelist = new edge_t[E];

    for (i=0; i<E; i++)
    {
      edgelist[i].v1 = e[(i*2)+0];
      edgelist[i].v2 = e[(i*2)+1];
      edgelist[i].w  = w[i];
    }
  }

  void sortEdgeList ()
  {
    unsigned int  i, t;

    fprintf (stdout, "Sorting edges...\0");
    t = clock ();

    std::sort<edge_t*> (&edgelist[0], &edgelist[E], MSTree::edgeComp);

    fprintf (stdout, "\rSorted edges in %0.2fs\n\0", (float)(clock () - t)/CLOCKS_PER_SEC);
  }

  void cleanEdgeList ()
  {
    if (edgelist != NULL)
    {
      delete[] edgelist;
      edgelist = NULL;
    }
  }

public:
  MSTree (
    Tree **tree,
    const unsigned int vertexCount,
    const unsigned int *edgeList,
    const float *edgeWeights,
    const unsigned int edgeCount,
    const bool verboseOutput)
    : V(vertexCount),
      E (edgeCount),
      e (edgeList),
      w (edgeWeights),
      edgelist (NULL),
      TE (vertexCount-1),
      TC (2*vertexCount-1),
      mst (tree),
      verbose (verboseOutput)
  {
    allocEdgeList ();
  }

  virtual ~MSTree ()
  {
    cleanEdgeList ();
  }

  /**
   * initialise the algorithm variables for processing
   */
  virtual void init () = 0;

  /**
   * find the next edge and add it to the tree
   * @return true if an edge was added, false otherwise
   */
  virtual bool next () = 0;
};

class MSTBoruvka : public MSTree
{
protected:
  unsigned int  *qmap;
  unsigned int  cmx;

public:
  MSTBoruvka (
    Tree **tree,
    const unsigned int vertexCount,
    const unsigned int *edgeList,
    const float *edgeWeights,
    const unsigned int edgeCount,
    const bool verboseOutput)
    : MSTree (tree, vertexCount, edgeList, edgeWeights, edgeCount, verboseOutput),
      qmap (NULL),
      cmx (0)
  {}

  ~MSTBoruvka ()
  {}

  void init ()
  {
    allocTree ();
  //  allocEdgeList ();
    initialComponents ();

#ifdef OP_QMAP
    unsigned int  i;

    qmap = new unsigned int[V];

    for (i=0; i<V; i++)
    {
      qmap[i] = i;
    }
#endif

    cmx = V;
  }

  bool next ()
  {
    unsigned int s,d,i,j;
    unsigned int *m;
  //  unsigned int icmx;
    unsigned int sr,dr,tr[2];
    unsigned int rcnt;
    float trw;

    Tree::Node *c1, *c2;
    Tree::Node *n;

    unsigned int oc[2];

    if ((rcnt = (*mst)->countRoots ()) == 1)
    {
    //  n = (*mst)->getNodeIdx (0);
    //  (*mst)->setRootID (n->getRoot ()->getCID ());
      (*mst)->setRoot ();

#ifdef OP_QMAP
      delete[] qmap;
#endif
      return false;
    }

    // map
    m = new unsigned int[TC];
    memset (m, 0, TC*sizeof (unsigned int));

  //  icmx = cmx;

    if (verbose)
    {
      fprintf (stdout, "\rmapping...\0");
    }

    // for each component
    for (i=0; i<(*mst)->length (); i++)
    {
      // check this component has not been mapped
      if (m[i])  {continue;}

      tr[0]=tr[1]=0;
      trw=FLT_MAX;

      // get the shortest edge connecting this component to another
      for (j=0; j<E; j++)
      {
        // get the src and dst of the edge
        s = e[(j<<1)+0];
        d = e[(j<<1)+1];

        // get the components of these vertices
#ifdef OP_QMAP
        sr = qmap[s];
        dr = qmap[d];
#else
        c1 = mst->findNode (s);
        c2 = mst->findNode (d);

      //  if ((c1==NULL)||(c2==NULL))
      //    continue;

      //  if ((c1->getRoot ()==NULL)||(c2->getRoot ()==NULL))
      //    continue;

        sr = c1->getRoot ()->getCID ();
        dr = c2->getRoot ()->getCID ();
#endif
        // if they are not our component, skip
        if ((sr!=i)&&(dr!=i))
        {
          continue;
        }

        // if the edge doesn't link components, skip
        if (sr==dr)
        {
          continue;
        }

    //    fprintf (stdout, "\rMST[Brv] %03i/%03i [%i/%i]...\0", cmx, mst->length (), j, E);

        // if it does, store it as the current minimum
        if (trw > w[j])
        {
          tr[0] = sr;
          tr[1] = dr;
          trw = w[j];
        }

        // early quit
        if (trw==0.0f)
        {
          break;
        }
      }

      // merge the components
      n = (*mst)->find (cmx);
      c1 = (*mst)->find (tr[0]);
      c2 = (*mst)->find (tr[1]);

      if ((n==NULL)||(c1==NULL)||(c2==NULL))
      {
        break;
      }

      n->setLeft (c1);
      n->setRight (c2);

      m[tr[0]] = 1;
      m[tr[1]] = 1;
    //  m[cmx] = 1;

      if (verbose)
      {
        fprintf (stdout, "\rCC: %i/%i [%i]\0", i, (*mst)->length (), (*mst)->countRoots ());
      }

#ifdef OP_QMAP
      oc[0] = tr[0];
      oc[1] = tr[1];

      for (j=0; j<V; j++)
      {
        if ((qmap[j]==oc[0])||(qmap[j]==oc[1]))
          qmap[j]=cmx;
      }
#endif
      if (++cmx > (*mst)->length ()-1)
      {
        break;
      }
    }

    if (verbose)
    {
      fprintf (stdout, "\n\0");
    }

    // one step complete
    delete[] m;

    // check if the tree is done at the beginning of the function
    return true;
  }
};

template<class T>
class NNTree : public MSTree
{
protected:
  unsigned int *qmap;
  unsigned int cmx;
  const T *verts;
  T *avgs;

  float (*w)(const T *v1, const T *v2);

public:
  NNTree (
    Tree **tree,
    const unsigned int vertexCount,
    const bool verboseOutput,
    float (*weightFunc)(const T *v1, const T *v2),
    const T *vertices)
    : MSTree (tree, vertexCount, NULL, NULL, 0, verboseOutput),
      qmap (NULL),
      cmx (0),
      w (weightFunc),
      verts (vertices)
  {}

  ~NNTree ()
  {}

  void init ()
  {
    unsigned int  i;

    allocTree ();
  //  allocEdgeList ();
    initialComponents ();

#ifdef OP_QMAP
    qmap = new unsigned int[V];

    for (i=0; i<V; i++)
    {
      qmap[i] = i;
    }
#endif

    cmx = V;

    avgs = new T[TC*5];
    memset (avgs, 0, TC*5*sizeof (T));

    for (i=0; i<V; i++)
    {
      avgs[(i*5)+0] = verts[(i*5)+0];
      avgs[(i*5)+1] = verts[(i*5)+1];
      avgs[(i*5)+2] = verts[(i*5)+2];
      avgs[(i*5)+3] = verts[(i*5)+3];
      avgs[(i*5)+4] = verts[(i*5)+4];
    }
  }

  bool next ()
  {
    unsigned int s,d,i,j;
  //  unsigned int *m;
    unsigned int sr,dr,tr[2];
    unsigned int rcnt;
    float trw, wt, mg1, mg2;

    Tree::Node *c1, *c2;
    Tree::Node *n;

    unsigned int oc[2];

    if ((rcnt = (*mst)->countRoots ()) == 1)
    {
    //  n = (*mst)->getNodeIdx (0);
    //  (*mst)->setRootID (n->getRoot ()->getCID ());
      (*mst)->setRoot ();

#ifdef OP_QMAP
      delete[] qmap;
#endif
      delete[] avgs;
      return false;
    }

    // map

    // for each component
    for (i=0; i<V; i++)
    {
      sr = qmap[i];

      tr[0]=tr[1]=0;
      trw=FLT_MAX;

      // find the nearest component to sr
      for (j=0; j<V; j++)
      {
        dr = qmap[j];

        // if they are in our component, skip
        if (dr==sr)
        {
          continue;
        }

        // store it as the current minimum
      //  wt = w (&verts[i*5], &verts[j*5]);
        wt = w (&avgs[sr*5], &avgs[dr*5]);

        if (trw > wt)
        {
          tr[0] = sr;
          tr[1] = dr;
          trw = wt;
        }

        // early quit
        if (trw==0.0f)
        {
          break;
        }
      }

      // merge the components
      n = (*mst)->find (cmx);
      n->setWeight (trw);
      c1 = (*mst)->find (tr[0]);
      c2 = (*mst)->find (tr[1]);

      if ((n==NULL)||(c1==NULL)||(c2==NULL))
      {
        break;
      }

      // set the left to be the lower (darker) of the two values
      for (j=0, mg1=0.0f, mg2=0.0f; j<5; j++)
      {
        mg1 += avgs[(tr[0]*5)+j]*avgs[(tr[0]*5)+j];
        mg2 += avgs[(tr[1]*5)+j]*avgs[(tr[1]*5)+j];
      }

      mg1 = sqrt (mg1);
      mg2 = sqrt (mg2);

      if (mg1 < mg2)
      {
        n->setLeft (c1);
        n->setRight (c2);
      }
      else
      {
        n->setRight (c1);
        n->setLeft (c2);
      }

    //  m[tr[0]] = 1;
    //  m[tr[1]] = 1;
    //  m[cmx] = 1;

      if ((verbose) && (i%100==0))
      {
        fprintf (stdout, "\rNNC: %i/%i [%i]\0", i, (*mst)->length (), (*mst)->countRoots ());
      }

#ifdef OP_QMAP
      oc[0] = tr[0];
      oc[1] = tr[1];

      for (j=0; j<V; j++)
      {
        if ((qmap[j]==oc[0])||(qmap[j]==oc[1]))
          qmap[j]=cmx;
      }
#endif
      // update the average colour
      avgs[(cmx*5)+0] = (avgs[(tr[0]*5)+0] + avgs[(tr[1]*5)+0])/2.0f;
      avgs[(cmx*5)+1] = (avgs[(tr[0]*5)+1] + avgs[(tr[1]*5)+1])/2.0f;
      avgs[(cmx*5)+2] = (avgs[(tr[0]*5)+2] + avgs[(tr[1]*5)+2])/2.0f;
      avgs[(cmx*5)+3] = (avgs[(tr[0]*5)+3] + avgs[(tr[1]*5)+3])/2.0f;
      avgs[(cmx*5)+4] = (avgs[(tr[0]*5)+4] + avgs[(tr[1]*5)+4])/2.0f;

      if (++cmx > (*mst)->length ()-1)
      {
        break;
      }
    }

    if (verbose)
    {
      fprintf (stdout, "\n\0");
    }

    // one step complete
  //  delete[] m;

    // check if the tree is done at the beginning of the function
    return true;
  }
};

class MSTreeX
{
public:
   /**
   * Build a minimum spanning tree of a graph using Prims method
   * @param V number of vertices in the graph (actually, the number of unique vertices referenced by the edgelist)
   * @param edges edgelist for the graph
   * @param weights weights of the edges used to build the tree
   * @param E number of edges in the edgelist
   * @param mstEdges pointer-to-pointer to recieve the edge list of the minimum spanning tree
   * @param mstWeights pointer-to-pointer to recieve the weights of the edges in the MST
   * @param mstEdgeCnt number of edges in the minimum spanning tree
   */
  static void buildMST_Prims (
    const unsigned int V,
    unsigned int *edges,
    float *weights,
    const unsigned int E,
    unsigned int **mstEdges,
    float **mstWeights,
    unsigned int& mstEdgeCnt)
  {
    const unsigned int  TE = V-1;

    // Implementation of Prims algorithm for an MST
    // FIXME: we must use Kruskals (I don't remember why i wrote this, assume an error!)
    bool *tree;      // indicates if a vertex is in the tree or not
    unsigned int *e;
    float *w;
    unsigned int i, v, n, k, c;
    float minw;

    // alloc the space
    tree = new bool[V];
    e = new unsigned int[TE*2];
    w = new float[TE];
    memset (e, 0, TE*2*sizeof (unsigned int));

    // init the list
    for (i=0; i<V; i++)  {tree[i] = false; w[i] = FLT_MAX;}

    // build the tree
    v = 0;    // first node
    n = 0;    // next node
    k = 0;    // count of edges in the mst

    while (1)
    {
      // mark v as being in the tree
      tree[v] = true;

      // check if there are any vertices to be added
      for (i=0; i<V && tree[i]; i++)
        ;;

      if (i==V)    // the whole list was true
        break;

      e[(k<<1)+0] = v;

      // find the minimum edge involving v
      for (i=0, n=0, minw=FLT_MAX;
           i<E;
           i++)
      {
        // check one of the nodes of this edge is v
        if ((edges[(i<<1)+0] != v) && (edges[(i<<1)+1] != v))
          continue;

        // one of the edge nodes is v, the other is the interesting one
        if (edges[(i<<1)+0] !=v)
          c = edges[(i<<1)+0];
        else
          c = edges[(i<<1)+1];

        // check the node is not in the tree already
        if (tree[c])
          continue;

        // check the weight
        if (weights[i] < w[k])
        {
          w[k] = weights[i];
          e[(k<<1)+1] = c;
        }
      }

      // set v = n
      v = e[(k<<1)+1];

      if ((++k) > TE)
        break;
    }

    // return the tree
    *mstEdges = e;
    *mstWeights = w;
    mstEdgeCnt = k;

    // cleanup
    delete[] tree;
    tree = NULL;
  }

  /**
   * Build a minimum spanning tree of a graph using Kruskals method
   * @param V number of vertices in the graph (actually, the number of unique vertices referenced by the edgelist)
   * @param edges edgelist for the graph
   * @param weights weights of the edges used to build the tree
   * @param E number of edges in the edgelist
   * @param mstEdges pointer-to-pointer to recieve the edge list of the minimum spanning tree
   * @param mstWeights pointer-to-pointer to recieve the weights of the edges in the MST
   * @param mstEdgeCnt number of edges in the minimum spanning tree
   * @param componentMap receives an array of size V describing the tree (i.e. v[i] = parent (v[i]))
   */
  static void buildMST_Kruskals (
    const unsigned int V,
    unsigned int *edges,
    float *weights,
    const unsigned int E,
    unsigned int **mstEdges,
    float **mstWeights,
    unsigned int& mstEdgeCnt,
    unsigned int **componentMap,
    unsigned int& componentMapCnt)
  {
    const unsigned int TE = V-1;
    unsigned int *e, *c, s,d,i,k;
    unsigned int cmx;
    unsigned int t;
    float *w;

#if OP_KRUSKAL_QMAP
    unsigned int    *qmap;
    unsigned int    oc[2], j;
#endif

    // alloc the space
    c = new unsigned int[2*V-1];    // component assignments
    w = new float[TE];        // weights
    e = new unsigned int[TE*2];    // edge list
    memset (e, 0, TE*2*sizeof (unsigned int));

    // set all weights to max
    for (i=0; i<TE; i++)
    {
      w[i] = FLT_MAX;
    }

    // each vertex is intially its own component
    for (i=0; i<2*V-1; i++)
    {
      c[i] = i;
    }

#if OP_KRUSKAL_QMAP
    qmap = new unsigned int[V];
    memset (qmap, 0, V*sizeof (unsigned int));

    for (i=0; i<V; i++)
    {
      qmap[i] = i;
    }
#endif

    // maximum is the end of the initial vertices
    cmx = V;

    // sort the edges of the complete graph
    GraphUtils::sortEdgeList (edges, E, weights);

    fprintf (stdout, "Building MST (Kruskal)...\0");
    t = clock ();

    // get the start point
    for (i=0, k=0; i<E; i++)
    {
    //  fprintf (stdout, "\rBuilding MST %i/%i\0", i, E);

      // get the src and dst of the edge
      s = edges[(i<<1)+0];
      d = edges[(i<<1)+1];

      // check if the edge links two components
#if OP_KRUSKAL_QMAP
      if (qmap[s] == qmap[d])
        continue;
#else
      if (PMT_getRoot (c, V, s) == PMT_getRoot (c, V, d))
        continue;
#endif

      // edge does link two components, add it to the mst
      e[(k<<1)+0] = s;
      e[(k<<1)+1] = d;
      w[k] = weights[i];

#if OP_KRUSKAL_QMAP
      // overwrite components
      oc[0] = qmap[s];
      oc[1] = qmap[d];

      for (j=0; j<V; j++)
      {
        if ((qmap[j] == oc[0]) || (qmap[j] == oc[1]))
          qmap[j] = cmx;
      }
#else
      // update the component map
      c[PMT_getRoot (c,s)] = cmx;
      c[PMT_getRoot (c,d)] = cmx;
#endif

      // increment the component count
      cmx++;

      // next edge pls
      if (++k > TE)
        break;
    }

    *mstEdges = e;
    *mstWeights = w;
    *componentMap = c;
    mstEdgeCnt = TE;
    componentMapCnt = cmx;

#if OP_KRUSKAL_QMAP
    fprintf (stdout, "\rConstructed MST (Kruskal-QMAP) in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
    delete[] qmap;
#else
    fprintf (stdout, "\rConstructed MST (Kruskal) in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
#endif
  }

  static Tree* buildMST_Boruvka (
    const unsigned int V,
    unsigned int *edges,
    float *weights,
    const unsigned int E)
  {
    const unsigned int TE = V-1;  // total edges in the mst
    const unsigned int TC = 2*V-1;  // total components in the graph
    unsigned int s,d,i,j;
    unsigned int *m;
    unsigned int cmx;
    unsigned int t,sr,dr,tr[2];
    float trw;

    Tree *mst;
    Tree::Node *c1, *c2;
    Tree::Node *n;

#ifdef OP_QMAP
    unsigned int *qmap, oc[2];
    qmap = new unsigned int[V];

    for (i=0; i<V; i++)
    {
      qmap[i] = i;
    }
#endif

    // initialise the tree
    mst = new Tree ();
    mst->setLength (TC);

    // map
    m = new unsigned int[TC];
    memset (m, 0, TC*sizeof (unsigned int));

    // each component is intially its own parent (a forest)
    for (i=0; i<TC; i++)
    {
      mst->get (i)->setCID (i);
    }

    fprintf (stdout, "MST[Brv]...\0");
    t = clock ();

    // while there is more than one component
    for (cmx=V; mst->countRoots () > 1;)
    {
      fprintf (stdout, "\rMST[Brv] %03i/%03i...\0", cmx, mst->length ());
      // merge the list of unmerged components
      memset (m, 0, TC*sizeof (unsigned int));

      // for each component
      for (i=0; i<mst->length (); i++)
      {
        // check this component has not been mapped
        if (m[i])  {continue;}

      //  fprintf (stdout, "\rMST[Brv] %03i/%03i...\0", cmx, mst->length ());

        tr[0]=tr[1]=0;
        trw=FLT_MAX;

        // get the shortest edge connecting this component to another
        for (j=0; j<E; j++)
        {
          // initial test of the weight, quick reject
          if (trw<weights[j])
            continue;

          // get the src and dst of the edge
          s = edges[(j<<1)+0];
          d = edges[(j<<1)+1];

          // get the components of these vertices
#ifdef OP_QMAP
          sr = qmap[s];
          dr = qmap[d];
#else
          c1 = mst->findNode (s);
          c2 = mst->findNode (d);

        //  if ((c1==NULL)||(c2==NULL))
        //    continue;

        //  if ((c1->getRoot ()==NULL)||(c2->getRoot ()==NULL))
        //    continue;

          sr = c1->getRoot ()->getCID ();
          dr = c2->getRoot ()->getCID ();
#endif
          // if they are not our component, skip
          if ((sr!=i)&&(dr!=i))
          {
            continue;
          }

          // if the edge doesn't link components, skip
          if (sr==dr)
          {
            continue;
          }

      //    fprintf (stdout, "\rMST[Brv] %03i/%03i [%i/%i]...\0", cmx, mst->length (), j, E);

          // if it does, store it as the current minimum
          if (trw > weights[j])
          {
            tr[0] = sr;
            tr[1] = dr;
            trw = weights[j];
          }

          // early quit
          if (trw==0.0f)
          {
            break;
          }
        }

        // merge the components
        n = mst->find (cmx);
        c1 = mst->find (tr[0]);
        c2 = mst->find (tr[1]);

        if ((n==NULL)||(c1==NULL)||(c2==NULL))
        {
          break;
        }

        n->setLeft (c1);
        n->setRight (c2);

        m[tr[0]] = 1;
        m[tr[1]] = 1;

#ifdef OP_QMAP
        oc[0] = tr[0];
        oc[1] = tr[1];

        for (j=0; j<V; j++)
        {
          if ((qmap[j]==oc[0])||(qmap[j]==oc[1]))
            qmap[j]=cmx;
        }
#endif

        if (++cmx > mst->length ()-1)
        {
          break;
        }
      }
    }

    fprintf (stdout, "\rMST[Brv] in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);

    fprintf (stdout, "Tree with %i roots\n\0", mst->countRoots ());
    // finish
    delete[] m;
#ifdef OP_QMAP
    delete[] qmap;
#endif

    if (mst->countRoots () == 1)
    {
      mst->setRoot ();
    //  n = mst->get (0);
    //  mst->setRootID (n->root ()->getCID ());
    }

    return mst;
  }
};
