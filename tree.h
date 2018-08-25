#include <map>

class Tree
{
public:
  typedef enum
  {
    PreOrder,
    InOrder,
    PostOrder
  } SearchOrder;

  class Node
  {
  protected:
    unsigned int cid; // component id
    float w; // weight
    Node *l, *r, *p;

  public:
    Node ()
      : w (0.0f), cid (0), l (NULL), r (NULL), p (NULL)
    {}

    ~Node ()
    {}

    void setWeight (const float newWeight)
    {
      w = newWeight;
    }

    void setCID (const unsigned int newCID)
    {
      cid = newCID;
    }

    const float weight () const
    {
      return w;
    }

    const unsigned int getCID () const
    {
      return cid;
    }

    Node *parent () const
    {
      return p;
    }

    Node *left () const
    {
      return l;
    }

    Node *right () const
    {
      return r;
    }

    void setLeft (Node *newLeft)
    {
      l = newLeft;

      if (l != NULL)
        l->p = this;
    }

    void setRight (Node *newRight)
    {
      r = newRight;

      if (r != NULL)
        r->p = this;
    }

    void setParent (Node *newParent)
    {
      p = newParent;
    }

    // misc
    const Node *getRoot () const
    {
      const Node *c = this;

      while (c->parent () != NULL)
      {
        c = c->parent ();
      }

      return c;
    }

    const bool isLeaf () const
    {
      return ((r == NULL) && (l == NULL));
    }

    /**
     * height (distance from leaf)
     */
    const unsigned int height () const
    {
      if ((l == NULL) && (r == NULL))
      {
        return 0;
      }
      else
      {
        unsigned int  ld, rd;

        ld = rd = 0;

        if (l != NULL)
          ld = left ()->height ();
        if (r != NULL)
          rd = right ()->height ();

        return 1 + ((ld<rd)?rd:ld); //std::max (ld, rd);
      }
    }

    const unsigned int depth () const
    {
      if (p == NULL)
      {
        return 0;
      }
      else
      {
        return 1 + p->depth ();
      }
    }

    const int dist (const Node* t) const  // distance to anscestor (t)
    {
      const Node* n;
      int l = 0;

      n = this;

      while (n != t)
      {
        if (n->parent () == NULL)
        {
          l = -1;
          break;
        }

        n = n->parent ();
        l++;
      }

      return l;
    }

    /**
     * @return number of nodes in this subtree
     * FIXME: return the number of leaves under this node
     */
    const unsigned int size () const
    {
      unsigned int  s;

      s = 1;

      if (left () != NULL)
      {
        s += left ()->size ();
      }
      if (right () != NULL)
      {
        s += right ()->size ();
      }

      return s;
    }

    /**
     * count the leaves under a node
     */
    const unsigned int countLeaves () const
    {
      unsigned int  s = 0;

      if (isLeaf ())
        return 1;

      if (left () != NULL)
        s += left ()->countLeaves ();
      if (right () != NULL)
        s += right ()->countLeaves ();

      return s;
    }

    /**
     * test if p is an ancestor of this node
     * @param p node to test
     * @return true, if p is an ancestor (appears on a path from this to the root)
     */
    const bool isAncestor (const Node *p) const
    {
      const Node *n = this;

      while (n != p)
      {
        n = n->parent ();

        if (n == NULL)
          return false;
      }

      return true;
    }
  };

protected:
  unsigned int len;
  Node *ns;
  const Node *r;

  // quick lookup table
  std::map<unsigned int, Node*>  lkp;

  void allocNodes ()
  {
    ns = new Node[len];
    memset (ns, 0, len*sizeof (Node));
  }

  void cleanNodes ()
  {
    if (ns != NULL)
    {
      delete[] ns;
      ns = NULL;
    }
    len = 0;
    r = NULL;
  }

  /**
   * iterates through the node list and replaces references
   * to p as parent with references to np
   */
  void setParentRefs (const Node *p, const Node *np)
  {
    unsigned int  i;

    for (i=0; i<length (); i++)
    {
      if (ns[i].parent () == p)
        ns[i].setParent (const_cast<Node*>(np));
    }
  }

public:
  Tree ()
    : ns (NULL), r (NULL), len (0)
  {}

  virtual ~Tree ()
  {
    cleanNodes ();
  }

  void setLength (const unsigned int newLength)
  {
    cleanNodes ();
    len = newLength;
    allocNodes ();
  }

  void setRoot ()
  {
    r = ns[0].getRoot ();
  }

  const unsigned int length () const
  {
    return len;
  }

  const Node *root () const
  {
    return r;
  }

  // search functions
  void buildLookupTable ()
  {
    unsigned int i;
    Node *n;

    for (i=0; i<length (); i++)
    {
      n = get (i);

      lkp[n->getCID ()] = n;
    }
  }

  Node *find (const unsigned int componentID)
  {
    unsigned int i;

    if (lkp.size () > 0)
    {
      return lkp[componentID];
    }
    else
    {
      for (i=0; i<length (); i++)
      {
        if (ns[i].getCID () == componentID)
          return &ns[i];
      }

      return NULL;
    }
  }

  Node *get (const unsigned int arrayIndex) const
  {
    return &ns[arrayIndex];
  }

  // factory function
  // lol wtf is PMT Rep?
  static Tree* fromPMTRep (const unsigned int *pmt, const unsigned int V)
  {
    const unsigned int BAD_VAL = 0xFFFFFFFF;
    unsigned int q, rq, lq, i, j;
    Tree *t;
    Node *n;

    t = new Tree ();
    t->setLength (V);

    q = GraphUtils::PMT_getRoot (pmt, pmt[1]);
  //  t->rid = q;

  //  or = new unsigned int[V];
  //  memset (or, 0xFF, V*sizeof (unsigned int));

  //  GraphUtils::PMT_getBreadthFirstOrdering (pmt, V, q, or);

    // add all the nodes
    for (i=0; i<V; i++)
    {
      n = t->get (i);
      n->setCID (i);
      n->setWeight (FLT_MAX);
    }

    // set the left and right nodes
    for (i=0; i<t->length (); i++)
    {
      // get the node for this component
    //  n = t->findNode (i);
      n = t->get (i);

      if (n == NULL)
        break;

      // get the left and right nodes
      for (j=0, rq=BAD_VAL, lq=BAD_VAL; j<V; j++)
      {
        if ((pmt[j] == n->getCID ()) && (pmt[j] != j))
        {
          if (lq == BAD_VAL)
          {
            lq = j;
          }
          else
          {
            rq = j;
            break;
          }
        }
      }

      if (lq != BAD_VAL)
      {
        n->setLeft (t->find (lq));
      }

      if (rq != BAD_VAL)
      {
        n->setRight (t->find (rq));
      }

    }

    // build the weights
  //  t->buildWeights ();
    // cleanup
  //  delete[] or;

    t->setRoot ();

    // return
    return t;
  }

  // utils
  /**
   * get an ordering of the tree
   * FIXME: if we passed the root to this function, it could be static and more versatile.
   * @param or stl list of node pointers, which contains (pointers to) the nodes in the order requested)
   * @param ordering, search order requested. NB: only InOrder ordering is supported for now.
   */
  void getDepthFirstOrdering (std::list<const Node*>& or, SearchOrder ordering)
  {
    std::list<const Node*>  s;
    const Node        *n;

    s.push_front (root ());

    switch (ordering)
    {
      case PreOrder:
      {
        break;
      }
      case InOrder:
      {
        while (!s.empty ())
        {
          n = s.front ();
          s.pop_front ();

          // process left
          if (n->left () != NULL)
          {
            // check if the left child is not in the ordering
            if (std::find (or.begin (), or.end (), n->left ()) == or.end ())
            {
              // no, push everything on and restart
              s.push_front (n);
              s.push_front (n->left ());
              continue;
            }
          }

        //  process the centre node
        //  {
            // add it to the list
      //    if (std::find (or.begin (), or.end (), n) == or.end ())
            or.push_back (n);
        //  }

          // right node
          if (n->right () != NULL)
          {
            s.push_front (n->right ());
          }
        }
        break;
      }
      case PostOrder:
      {
        break;
      }
    }
  }

  /**
   * basic node weighting function
   * sets node weights to the number of nodes in the subtree with that node at the root.
   */
  void buildWeights ()
  {
    unsigned int  i;
    Node      *n;

    for (i=0; i<length (); i++)
    {
      n = get (i);

      n->setWeight ((float)n->size ());
    }
  }

  /**
   * count the number of nodes with no parent in the array
   * @return number of root nodes encountered
   */
  const unsigned int countRoots () const
  {
    unsigned int  i, c;

    for (i=0,c=0; i<length (); i++)
    {
      if (ns[i].parent () == NULL)
        c++;
    }

    return c;
  }

  /**
   * count the number of leaves in the array.
   * NB: this iterates over the array, so store the result
   * instead of repeated calls.
   * @return the number of nodes in the array for which isLeaf () returns true
   */
  const unsigned int countLeaves () const
  {
    unsigned int  i, c;

    for (i=0, c=0; i<length (); i++)
    {
      if (ns[i].isLeaf ())
        c++;
    }

    return c;
  }

  /**
   * count the number of nodes referencing n as a parent
   * @param n node to check
   * @return the number of nodes referencing n as a parent
   */
  const unsigned int countParentRefs (const Node *n) const
  {
    unsigned int  i, cnt;

    for (i=0, cnt=0; i<length (); i++)
    {
      if (ns[i].parent () == n)
        cnt++;
    }

    return cnt;
  }

  /**
   * returns true if this node is a parent to an existing node
   * @param n node to test
   */
  const bool isParent (const Node *n) const
  {
    unsigned int  i;

    for (i=0; i<length (); i++)
    {
      if (ns[i].parent () == n)
        return true;
    }

    return false;
  }

  /**
   * get a list of the internal leaves
   */
  void getInternalLeafList (std::list<const Node*>& ls)
  {
    unsigned int  i;
    const Node    *n;

    ls.clear ();

    for (i=0; i<length (); i++)
    {
      n = &ns[i];

      if ((n->isLeaf ()) && (isParent (n)))
      {
        ls.push_back (n);
      }
    }
  }

  /**
   * write a tree to a filestream
   * @param f filestream to write to
   */
  void write (const char *cacheDir, const char *fname)
  {
    unsigned int cid[4];
    unsigned int i, l, t;
    float w;
    Node *n;
    FILE *f;
    char path[MAX_PATH];

    memset (path, 0, MAX_PATH);
    sprintf (path, "%s%s.tree\0", cacheDir, fname);

    fprintf (stdout, "Writing tree...\0");
    t = clock ();

    f = fopen (path, "wb\0");

    if (f == NULL)
      return;

    // write the header
#define WRITE_MINIMAL_TREE
#ifdef WRITE_MINIMAL_TREE
    l = root ()->size ();
#else
    l = len;
#endif
    fwrite (&l, sizeof (unsigned int), 1, f);

    // write the root id
    cid[0] = root ()->getCID ();
    fwrite (&cid[0], sizeof (unsigned int), 1, f);

    // write the nodes
    for (i=0; i<length (); i++)
    {
      n = get (i);
#ifdef WRITE_MINIMAL_TREE
      if (!n->isAncestor (root ()))
        continue;
#endif

      cid[0] = n->getCID ();
      if (n->left () == NULL)   {cid[1] = 0xFFFFFFFF;} else {cid[1] = n->left ()->getCID ();}
      if (n->right () == NULL)  {cid[2] = 0xFFFFFFFF;} else {cid[2] = n->right ()->getCID ();}
      if (n->parent () == NULL) {cid[3] = 0xFFFFFFFF;} else {cid[3] = n->parent ()->getCID ();}
      w = n->weight ();

      fwrite (cid, sizeof (unsigned int), 4, f);
      fwrite (&w, sizeof (float), 1, f);
    }

    // close
    fclose (f);

    fprintf (stdout, "\rWrote tree in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
  }

  /**
   * read a file from a filestream
   * @param f filestream to read from
   */
  static Tree* read (const char *cacheDir, const char *fname)
  {
    unsigned int l, rcid, cid[4];
    unsigned int i, *cids;
    unsigned int tm;
    float w;
    Tree *t;
    Node *n;
    FILE *f;
    char path[MAX_PATH];
#define OP_NODE_READ_MAP
#ifdef OP_NODE_READ_MAP
    std::map<unsigned int, Node*> lkp;  // lookup accl
#endif
    fprintf (stdout, "Reading tree...\0");
    tm = clock ();

    memset (path, 0, MAX_PATH);
    sprintf (path, "%s%s.tree\0", cacheDir, fname);

    f = fopen (path, "rb\0");

    if (f == NULL)
    {
      fprintf (stdout, "\rFile not found '%s'\n\0", fname);
      return NULL;
    }

    // get the length
    fread (&l, sizeof (unsigned int), 1, f);
    // get the root id
    fread (&rcid, sizeof (unsigned int), 1, f);

    fprintf (stdout, "\rReading tree [%i]...\0", l);

    t = new Tree ();
    t->setLength (l);

    // store the read cid of the left, right and parent for each node
    cids = new unsigned int[t->length ()*3];
    memset (cids, 0xFF, t->length ()*3*sizeof (unsigned int));
#ifdef OP_NODE_READ_MAP
  //  np = new Node*[t->length ()];
  //  memset (np, 0, t->length ()*sizeof (Node*));
#endif
    // read nodes
    for (i=0; i<t->length (); i++)
    {
      n = t->get (i);

      fread (cid, sizeof (unsigned int), 4, f);
      fread (&w, sizeof (float), 1, f);

      // set the cid and weight
      n->setCID (cid[0]);
      n->setWeight (w);

      // copy the lrp cids
      cids[(i*3)+0] = cid[1];
      cids[(i*3)+1] = cid[2];
      cids[(i*3)+2] = cid[3];

#ifdef OP_NODE_READ_MAP
    //  if (cid[0] > t->length ())
    //    fprintf (stdout, "ERR: CID > length () [%i > %i]\n\0", cid[0], t->length ());
    //  else
    //    np[cid[0]] = n;
      lkp[cid[0]] = n;
#endif
    }

    // close the file
    fclose (f);

    // once all the data is read we can set the pointers correctly
    for (i=0; i<t->length (); i++)
    {
      n = t->get (i);

      n->setLeft (NULL);
      n->setRight (NULL);
      n->setParent (NULL);

#ifdef OP_NODE_READ_MAP
      // fast lookup
      if (cids[(i*3)+0] != 0xFFFFFFFF)  {n->setLeft (lkp[cids[(i*3)+0]]);}
      if (cids[(i*3)+1] != 0xFFFFFFFF)  {n->setRight (lkp[cids[(i*3)+1]]);}
      if (cids[(i*3)+2] != 0xFFFFFFFF)  {n->setParent (lkp[cids[(i*3)+2]]);}
#else
      // slow lookup
      if (cids[(i*3)+0] != 0xFFFFFFFF)  {n->setLeft (t->find (cids[(i*3)+0]));}
      if (cids[(i*3)+1] != 0xFFFFFFFF)  {n->setRight (t->find (cids[(i*3)+1]));}
      if (cids[(i*3)+2] != 0xFFFFFFFF)  {n->setParent (t->find (cids[(i*3)+2]));}
#endif
    }

    // cleanup
    delete[] cids;
#ifdef OP_NODE_READ_MAP
  //  delete[] np;
#endif

    // set the root
    t->r = t->find (rcid);

    fprintf (stdout, "\rRead tree [%i] in %0.2fs\n\0", t->length (), (float)(clock ()-tm)/CLOCKS_PER_SEC);

    // finished
    return t;
  }

};

class TreePruner
{
protected:
  const bool verboseOutput;
  Tree *tree;

  void removeNodeRefs (Tree::Node *n)
  {
    Tree::Node  *p, *l, *r;

    p = n->parent ();
    l = n->left ();
    r = n->right ();

    // this would be a major fuck up
    if (p == NULL)  {return;}

    // set our childrens parent as our parent
  //  if (l != NULL)  {l->setParent (p);}
  //  if (r != NULL)  {r->setParent (p);}
  //  setParentRefs (n, p);

    // find which node we were in the parent
    if (p->left () == n)  {p->setLeft (NULL);}
    if (p->right () == n)  {p->setRight (NULL);}

    // remove the parent ref from this node to drop it from the tree
    n->setParent (NULL);
  }

  void updateLCMap (DataComponentMap *dcm, const unsigned int ocid, const unsigned int ncid)
  {
    unsigned int  i;
    Tree::Node    *p, *n;

    if (dcm != NULL)
    {
      p = tree->find (ocid);
    //  p = lkp[ocid];

      // replace the ids
      for (i=0; i<dcm->length (); i++)
      {
        n = tree->find (dcm->getDataComponent (i));
      //  n = lkp[cmpMap[i]];

      //  if ((cmpMap[i] == ocid) || (n->isAncestor (p)))
      //    cmpMap[i] = ncid;
        if ((dcm->getDataComponent (i) == ocid) || (n->isAncestor (p)))
          dcm->setDataComponent (i, ncid);
      }
    }
  }

public:
  TreePruner (Tree *t, const bool verbose)
    : tree (t), verboseOutput (verbose)
  {}

  virtual ~TreePruner ()
  {}

  virtual void init ()
  {
    tree->buildLookupTable ();
  }

  virtual void prune (const unsigned int minVolume, DataComponentMap *dcm) = NULL;

  const bool verbose ()
  {
    return verboseOutput;
  }
};

class TreePruner_Linear : public TreePruner
{
public:
  TreePruner_Linear (Tree *tree, const bool verbose)
    : TreePruner (tree, verbose)
  {}

  ~TreePruner_Linear ()
  {}

  void init ()
  {
    TreePruner::init ();
  }

  /**
   * Remove internal nodes in the tree which have less than some threshold value.
   * Leaf nodes are preserved; the parent pointer for leaf nodes connects to the
   * appropriate preserved internal node, while the left and right child of the
   * p-node remains NULL. This allows upward traversal from the leaves, while
   * reducing the overall size of the tree.
   * @param minVolume minimum volume of the node to be preserved.
   */
  void prune (const unsigned int minVolume, DataComponentMap *dcm)
  {
    unsigned int i, cnt, t;
    Tree::Node *n, *p;
    unsigned int *ss;

    fprintf (stdout, "Pruning nodes 0 < size () < %i\0", minVolume);
    t = clock ();

    // cache the sizes of the nodes (this data is not constant when we modify the tree)
    ss = new unsigned int[tree->length ()];
    memset (ss, 0, tree->length () * sizeof (unsigned int));

    for (i=0; i<tree->length (); i++)
    {
    //  ss[i] = ns[i].size ();
      ss[i] = tree->get (i)->countLeaves ();
    }

    for (i=0, cnt=0; i<tree->length (); i++)
    {
      n = tree->get (i);//&ns[i];

      if (verbose ())
      {
        fprintf (stdout, "\rPruning %i/%i [%i]\0", i, tree->length (), cnt);
      }

      if (n->parent () == NULL)  // don't remove roots
        continue;

      if (/*(!n->isLeaf ()) &&*/ (ss[i] < minVolume))//(n->size () < minVolume))
      {
        // store the parent
        p = n->parent ();

        // remove this node from the tree
        removeNodeRefs (n);

        // update the component id map, if we have one
        updateLCMap (dcm, n->getCID (), p->getCID ());

        // next
        cnt++;
      }
    }

    fprintf (stdout, "\rPruned %i nodes with 0 < size () < %i in %0.2fs\n\0", cnt, minVolume, (float)(clock ()-t)/CLOCKS_PER_SEC);

    delete[] ss;

    // FIXME: realloc the space and rebuild the tree with smaller memory footprint
  }
};

class TreePruner_DFS : public TreePruner
{
protected:
  DataComponentMap  *dcm;
  unsigned int    mv;

  /**
   * perform a depth first search on the tree and remove trees with less than minVolume
   * children. This should be a pre-order dfs, so we test the node before processing the
   * childre, otherwise we will remove all the nodes.
   * @param r node to test
   */
  void pruneNode (Tree::Node *r)
  {
    // test this node
  //  if (r->size () < mv)
  //  if (r->countLeaves () < mv)
    if (r->weight () < 0.01f)
    {
      unsigned int  ocid;

      if (verbose ())
      {
        fprintf (stdout, "\rPruning node %i\0", r->getCID ());
      }

      // rather than removing this node, remove its children,
      // this prevents a node having its left child removed, but its
      // right child remaining. i.e. we lose a leaf.

      if (r->left () != NULL)
      {
        ocid = r->left ()->getCID ();
        removeNodeRefs (r->left ());
        updateLCMap (dcm, ocid, r->getCID ());
      }

      if (r->right () != NULL)
      {
        ocid = r->right ()->getCID ();
        removeNodeRefs (r->right ());
        updateLCMap (dcm, ocid, r->getCID ());
      }
      return;
    }

    // else operate on the children
    if (r->left () != NULL)    {pruneNode (r->left ());}
    if (r->right () != NULL)  {pruneNode (r->right ());}
  }

public:
  TreePruner_DFS (Tree *tree, const bool verbose)
    : TreePruner (tree, verbose)
  {}

  ~TreePruner_DFS ()
  {}

  void init ()
  {
    TreePruner::init ();
  }

  void prune (const unsigned int minVolume, DataComponentMap *dcmap)
  {
    const Tree::Node  *r;

    mv = minVolume;
    dcm = dcmap;

    if (verbose ())
      fprintf (stdout, "\nPruning node...\0");

    // never operate on the root
    r = tree->root ();

    if (r->left () != NULL)   {pruneNode (r->left ());}
    if (r->right () != NULL)  {pruneNode (r->right ());}

  }
};

