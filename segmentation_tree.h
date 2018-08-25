#include <utility>
#include <map>

#define CC_IMG_DIR  "ccimgs"

class TreePlotter
{
public:
  static const unsigned int VALID = 0x00000001; ///< defunct
  static const unsigned int DRAW  = 0x00000002; ///< whether a node (and its children) should be drawn
  static const unsigned int CLPS  = 0x00000004; ///< whether a node is collapsed

  typedef struct DParm_s
  {
    unsigned int state;
    float pt[2];
    float w; ///< width of the parent space occupied by this node
    float xr[2]; ///< sizeof the horizontal region occupied by this node
  } DParm;

public:
  typedef enum PlotMethod
  {
    DFS,
    INDENT,
    TR,
    GRUNDO,
    AMDAM
  };

  typedef std::map<unsigned int, DParm> TreePlot;
  typedef TreePlot::iterator TreePlotIterator;
  typedef float Point[2];

  static void getTreeCoords_Basic (Tree *tree, TreePlot *coords)
  {
    std::list<const Tree::Node*>::iterator it;
    std::list<const Tree::Node*> or;
    unsigned int p, i;
    float y;
    const Tree::Node *n;

    tree->getDepthFirstOrdering (or, Tree::InOrder);

    // build the width of each subtree at each level
    for (it=or.begin (), p=0; it!=or.end (); it++)
    {
      n = *it;
      y = (float)(n->depth ());

    //  if ((*coords).find ((*it)->getCID ()) != (*coords).end ())
    //  {
    //    i = (*coords)[(*it)->getCID ()].state;
    //  }
    //  else
      {
        i = VALID|DRAW;
      }

    //  (*coords)[(*it)->getCID ()] = Point ((float)p, (float)((*it)->depth ()));
    //  (*coords)[(*it)->getCID ()] = Point ((float)p, (float)((*it)->weight ()));
    //  (*coords)[(*it)->getCID ()] = Point ((float)p, y);
      (*coords)[(*it)->getCID ()].state = i;
      (*coords)[(*it)->getCID ()].pt[0] = (float)p;
      (*coords)[(*it)->getCID ()].pt[1] = y;

      p++;
    }
  }
  static void getTreeCoords_Indented (Tree *tree, TreePlot *coords)
  {
  //  for (i=0; i<tree->getNodeIdx (i); i++)
    {

    }
  }
  static void getTreeCoords_TR (Tree *tree, TreePlot *coords)
  {
#if 0
    unsigned int i;
    const float sep = 1.0f;
    Tree::Node *n;      // parent of the nodes we are operating on
    Tree::Node *l, *lr;  // left subtree (right child iterates)
    Tree::Node *r, *rl;  // right subtree (left child iterates)
    float x;
    TreePlot tmp;

    // 1st: compute the offset of each node, left or right from its parent
    for (i=0; i<tree->length (); i++)
    {
      n = tree->get (i);  // get the node
    //  pt.second = n->depth ();      // get the level

      // enumerate its children
      l = n->left ();
    //  x = +sep;
      x = 0.0f;

      if (l != NULL)
      {
        for (lr=l; (lr=lr->right ()) != NULL; x-=sep)
          ;;

        tmp[l->getCID ()].first = x;
      }

      // do the right
      r = n->right ();
    //  x = -sep;
      x = 0.0f;

      if (r != NULL)
      {
        for (rl=r; (rl = rl->left ()) != NULL; x+=sep)
          ;;

        tmp[r->getCID ()].first = x;
      }
    }

    // 2nd: accumulate the offsets from the root to each node to get its actual position
    for (i=0; i<tree->length (); i++)
    {
      n = tree->get (i);

      // only operate on the x
      if (n->parent () != NULL)
      {
        if (n->parent ()->left () == n)
          x = -sep;
        else
          x = sep;
      }
      else
        x = 0.0f;

      // go to the root, mutating x as we travel (NB: r means root here)
      for (r=n; (r=r->parent ()) != NULL;)
      {
        x += tmp[r->getCID ()].first;
      }

      (*coords)[n->getCID ()] = Point ((float)x, (float)n->depth ());
    }
#endif
  }

  static void setNodeCoords_Grundo (Tree *tree, TreePlot *coords, const Tree::Node *n)
  {
    unsigned int p, lp, rp;

    if (n->isLeaf ())
      return;

    // post-order depth first processing
    if (n->left () != NULL)
      setNodeCoords_Grundo (tree, coords, n->left ());

    if (n->right () != NULL)
      setNodeCoords_Grundo (tree, coords, n->right ());

    // get the position of the left and right child
    // and position this node between them
    if (n->left () != NULL)
    {
      lp = (unsigned int)(*coords)[n->left ()->getCID ()].pt[0];
    }
    else
    {
      lp = 0xFFFFFFFF;
    }
    if (n->right () != NULL)
    {
      rp = (unsigned int)(*coords)[n->right ()->getCID ()].pt[0];
    }
    else
    {
      rp = 0xFFFFFFFF;
    }

    if ((lp != 0xFFFFFFFF) && (rp != 0xFFFFFFFF))
    {
      p = (lp + rp)/2;
    }
    else
    {
      if (lp != 0xFFFFFFFF)
        p = lp;
      else
        p = rp;
    }

    // update the plot
    (*coords)[n->getCID ()].pt[0] = (float)p;
  }

  static void getTreeCoords_Grundo (Tree *tree, TreePlot *coords)
  {
    // get the basic ordering
    getTreeCoords_Basic (tree, coords);

    // reorder the nodes based on child location
    setNodeCoords_Grundo (tree, coords, tree->root ());
  }

  static void setNodeCoords_Amsterdam (Tree *tree, TreePlot *coords, const Tree::Node *n, const unsigned int lx, const unsigned int rx)
  {
    const Tree::Node *l, *r;
    unsigned int ns, ls, rs;
    unsigned int lrng[2], rrng[2];

    // set the range of this node
    (*coords)[n->getCID ()].xr[0] = (float)lx;
    (*coords)[n->getCID ()].xr[1] = (float)rx;
    (*coords)[n->getCID ()].pt[0] = (float)lx + ((float)(rx-lx)/2.0f);

    // ignore leaves
    if (n->isLeaf ())
      return;

    l = n->left ();
    r = n->right ();

    if ((l == NULL) || (r == NULL))
    {
      if (l == NULL)  {(*coords)[r->getCID ()].w = 1.0f;}
      if (r == NULL)  {(*coords)[l->getCID ()].w = 1.0f;}
    }
    else // get the ratio of the sizes
    {
      // get the sizes of the nodes
      ns = n->size ()-1;
      ls = l->size ();
      rs = r->size ();

      (*coords)[l->getCID ()].w = (float)ls/(float)ns;
      (*coords)[r->getCID ()].w = (float)rs/(float)ns;
    }

    // process the children
    if (l != NULL)
    {
      lrng[0] = lx;
      lrng[1] = lx + ((rx-lx)*(*coords)[l->getCID ()].w);
      setNodeCoords_Amsterdam (tree, coords, l, lrng[0], lrng[1]);
    }

    if (r != NULL)
    {
      rrng[0] = rx - ((rx-lx)*(*coords)[r->getCID ()].w);
      rrng[1] = rx;
      setNodeCoords_Amsterdam (tree, coords, r, rrng[0], rrng[1]);
    }
  }

  static void getTreeCoords_Amsterdam (Tree *tree, TreePlot *coords)
  {
    TreePlotIterator  it;
    unsigned int    rng[2];
    unsigned int    x;

    // get the initial basic ordering
    getTreeCoords_Grundo (tree, coords);

    // get the range
    rng[0] = 0xFFFFFFFF;
    rng[1] = 0;

    for (it=coords->begin (); it!=coords->end (); it++)
    {
      x = (unsigned int)it->second.pt[0];

      if (rng[0]>x)  {rng[0]=x;}
      if (rng[1]<x)  {rng[1]=x;}
    }

    // set the width component of the nodes
    setNodeCoords_Amsterdam (tree, coords, tree->root (), rng[0], rng[1]);
  }

  static void getTreeCoords (Tree *tree, TreePlot *coords, PlotMethod pm)
  {
    switch (pm)
    {
    case DFS: getTreeCoords_Basic (tree, coords);    break;
    case INDENT: getTreeCoords_Indented (tree, coords);  break;
    case TR: getTreeCoords_TR (tree, coords);    break;
    case GRUNDO: getTreeCoords_Grundo (tree, coords);  break;
    case AMDAM: getTreeCoords_Amsterdam (tree, coords);  break;
    }
  }
};

class SegmentationTree
{
public:
  typedef enum LayoutStyle
  {
    TRAD,
    RADIAL
  };
  typedef enum DrawStyle
  {
    NODELINK,
    BLOCK,
    BLOCK_TEXTURES
  };

protected:
  const unsigned int  refTexture; ///< reference image

  // MST hierarchy
  Tree *hTree;
  DataComponentMap *dcm;
  TreePlotter::TreePlot *hPlot;
  float tp_rng[4];

  // component set
  ComponentSet *componentSet; // for images the data space is two dimensional

  // interactions
  unsigned int selectedComponent; ///< currently selected component
  Tree::Node *selectedNode; ///< currently selected node
  LayoutStyle layoutStyle;
  DrawStyle drawStyle;
  float split; ///< tree separation

  void createStipplePattern ()
  {
#if 0
    unsigned char  d[32*(32/8)];
    unsigned char  p;

    for(int i=0;i<32*(32/8);++i)
    {
      p = 0;

      for (int j=0; j<8; j++)
      {
        if (((float)(rand ())/RAND_MAX) > 0.5f)
          p |= 1<<j;
      }
      d[i]=p;

      srand (rand ());
    }
#else
    unsigned char  d[32*(32/8)];

    for (int y=0;y<32;y++)
    {
      for (int x=0;x<32/8;x++)
      {
        d[(y*32/8)+x] = 0;

        if (y%2 == 0)
        {
          d[(y*32/8)+x] = 0xCC;
        }
        else
        {
          d[(y*32/8)+x] = 0x33;
        }

      }
    }
#endif

    glPolygonStipple (d);
  }

  // cartesian to radial
  static void c2r (float *p)
  {
    const float  r = p[1];
    const float  th = p[0]*2.0f*(float)M_PI;

    p[0] = r * (cos (th) - sin (th));
    p[1] = r * (cos (th) + sin (th));
  }

  void drawCircleTree ()
  {
    const unsigned int  max_d = hTree->root ()->height ()+1;
    const unsigned int  max_j = 180;
    unsigned int    i, j;
    float        x, y;
    float  r, th;

    // draw circles at radius 0 to max_d
    for (i=1; i<max_d+1; i++)
    {
      r = (float)i/(float)(max_d);

    //  y = 1.0f-y;

      glBegin (GL_LINE_LOOP);
      glColor3f (0.8f, 0.8f, 0.8f);

      for (j=0; j<max_j; j++)
      {
        x = (float)j/(float)(max_j-1);

        // make it radial
        th = (float)(x*2.0*M_PI); // 2.0 * ;

        // rotate the point by the texture coord
        y = (r) * (cos (th) - sin (th));
        x = (r) * (cos (th) + sin (th));

        glVertex2f (x, y);
      }

      glEnd ();
    }
  }

  void setColour (const ComponentSet::Component *c, const ComponentSet::Component *t)
  {
    if (t == NULL)
    {
    //  glColor3fv (c->getAverage ());
      glColor3f (0,0,0);
    }
    else
    {
      unsigned int  i;
      float      r;
      const float    *cav, *tav;

      cav = c->getAverage ();
      tav = t->getAverage ();

      for (i=0, r=0.0f; i<4; i++)
      {
        r += (cav[i]-tav[i])*(cav[i]-tav[i]);
      }

      r = sqrt (r);

      if (r > 1.0f)
        r = 1.0f;

      r = 1.0f-r;

      r = pow (r, 2.0f);

      glColor3f (r, 0.0f, 0.0f);//, cav[1], cav[2]);
    }
  }

  void drawMSTNode (const Tree::Node *n, const bool selected)
  {
    ComponentSet::Component *cc, *sc;
    unsigned int cid;
    const float y_spc = 1.0f/(float)(hTree->root ()->height ()+1);//0.05f;//1.0f/(tp_rng[3]-tp_rng[1]);
    const float x_spc = 0.001f;//1.0f/(float)(hTree->length ()+2);
    Tree::Node *l, *r;
    float pt[2], pt2[4];
    const float *sav;

    cid = n->getCID ();
    cc = componentSet->find (cid);
    l = n->left ();
    r = n->right ();

    // set the selected state
    if (selectedNode == n)
    {
      glDisable (GL_POLYGON_STIPPLE);
      (bool)(selected) = true;
    }

    // set the stipple pattern for drawing unselected nodes
    if (selected)
    {
      glDisable (GL_POLYGON_STIPPLE);
    }
    else
    {
      glEnable (GL_POLYGON_STIPPLE);
    }

    if (selectedNode != NULL)
    {
      sc = componentSet->find (selectedNode->getCID ());
      sav = sc->getAverage ();
    }
    else
    {
      sc = NULL;
      sav = NULL;
    }

    pt[0] = (*hPlot)[cid].pt[0];  pt[1] = (*hPlot)[cid].pt[1];
    pt2[0] = (*hPlot)[cid].xr[0];  pt2[1] = pt[1];
    pt2[2] = (*hPlot)[cid].xr[1];  pt2[3] = pt[1]+y_spc;

    // draw the node
    if ((drawStyle == BLOCK) || (drawStyle == BLOCK_TEXTURES))
    {
      if (layoutStyle == RADIAL)
      {
        // block radial is a bit of work
        float  x, ptr[4];

        glColor3fv (cc->getAverage ());
      //  setColour (cc, sc);

        if ((drawStyle == BLOCK_TEXTURES) && (cc->hasTab ()))
        {
          glEnable (GL_TEXTURE_2D);
          glBindTexture (GL_TEXTURE_2D, cc->getTabID ());
          glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        }

        glBegin (GL_QUAD_STRIP);
        for (x=pt2[0]; x<pt2[2]; x+=x_spc)
        {
          ptr[0]=x;  ptr[1]=pt[1];
          ptr[2]=x;  ptr[3]=pt[1]+y_spc;

          c2r (ptr);    glTexCoord2f ((x-pt2[0])/(pt2[2]-pt2[0]), 0.0f);  glVertex2fv (&ptr[0]);
          c2r (&ptr[2]);  glTexCoord2f ((x-pt2[0])/(pt2[2]-pt2[0]), 1.0f);  glVertex2fv (&ptr[2]);
        }
        glEnd ();

        if (glIsEnabled (GL_TEXTURE_2D) == GL_TRUE)//((cc->hasTab ())
        {
          glDisable (GL_TEXTURE_2D);
        }
#define NODE_OUTLINES
#ifdef NODE_OUTLINES
        if (selected)//(isSelected (cid))
          glColor3f (0,0,0);
        else
          setColour (cc, sc);
        {
          glLineWidth (2);
      //    glColor3f (0,0,0);

          glBegin (GL_LINE_LOOP);
          ptr[0]=pt2[0];  ptr[1]=pt[1];
          ptr[2]=pt2[0];  ptr[3]=pt[1]+y_spc;
          c2r (ptr);        c2r (&ptr[2]);

          glVertex2fv (&ptr[0]);  glVertex2fv (&ptr[2]);

          for (x=pt2[0]; x<pt2[2]; x+=x_spc)
          {
            ptr[0]=x;  ptr[1]=pt[1]+y_spc;  c2r (ptr);
            glVertex2fv (&ptr[0]);
          }

          ptr[0]=pt2[2];  ptr[1]=pt[1];
          ptr[2]=pt2[2];  ptr[3]=pt[1]+y_spc;
          c2r (ptr);    c2r (&ptr[2]);

          glVertex2fv (&ptr[2]);  glVertex2fv (&ptr[0]);

          for (x=pt2[2]; x>pt2[0]; x-=x_spc)
          {
            ptr[0]=x;  ptr[1]=pt[1];  c2r (ptr);
            glVertex2fv (&ptr[0]);
          }

          glEnd ();
          glLineWidth (1);
        }
#endif
      }
      else
      {
        if ((drawStyle == BLOCK_TEXTURES) && (cc->hasTab ()))
        {
          glEnable (GL_TEXTURE_2D);
          glBindTexture (GL_TEXTURE_2D, cc->getTabID ());
          glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        }
        else
        {
          glColor3fv (cc->getAverage ());
        }

        glBegin (GL_QUADS);
        glTexCoord2f (0, 0);  glVertex2f (pt2[0]-x_spc, pt2[1]);
        glTexCoord2f (1, 0);  glVertex2f (pt2[2]-x_spc, pt2[1]);
        glTexCoord2f (1, 1);  glVertex2f (pt2[2]-x_spc, pt2[3]);
        glTexCoord2f (0, 1);  glVertex2f (pt2[0]-x_spc, pt2[3]);
        glEnd ();

        if (glIsEnabled (GL_TEXTURE_2D) == GL_TRUE)
        {
          glDisable (GL_TEXTURE_2D);
        }
#ifdef NODE_OUTLINES
        if (selected)//(isSelected (cid))
          glColor3f (0,0,0);
        else
          setColour (cc, sc);
        {
          glLineWidth (2);
        //  glColor3f (0,0,0);
          glBegin (GL_LINE_LOOP);
          glVertex2f (pt2[0]-x_spc, pt2[1]);  glVertex2f (pt2[2]-x_spc, pt2[1]);
          glVertex2f (pt2[2]-x_spc, pt2[3]);  glVertex2f (pt2[0]-x_spc, pt2[3]);
          glEnd ();
          glLineWidth (1);
        }
#endif
      }
    }
    else if (drawStyle == NODELINK)
    {
      if (layoutStyle == RADIAL)  {c2r (pt);}

    //  glColor3fv (cc->getAverage ());
    //  glColor3fv (cc->getAverage ());
      setColour (cc, sc);
      glBegin (GL_LINES);

      // draw a connecting edge to children
      if (l != NULL)
      {
        pt2[0] = (*hPlot)[l->getCID ()].pt[0];
        pt2[1] = (*hPlot)[l->getCID ()].pt[1];
        if (layoutStyle == RADIAL)  {c2r (pt2);}
        glVertex2f (pt[0], pt[1]);
        glVertex2f (pt2[0], pt2[1]);
      }

      if (r != NULL)
      {
        pt2[0] = (*hPlot)[r->getCID ()].pt[0];
        pt2[1] = (*hPlot)[r->getCID ()].pt[1];
        if (layoutStyle == RADIAL)  {c2r (pt2);}
        glVertex2f (pt[0], pt[1]);
        glVertex2f (pt2[0], pt2[1]);
      }

      glEnd ();

      // draw the node
    //  if (isSelected (cid))
      setColour (cc, sc);
      {
        glPointSize (12);
      //  glColor3f (1,1,1);
      //  glColor3f (0,0,0);
        glBegin (GL_POINTS);
        glVertex2f (pt[0], pt[1]);
        glEnd ();
        glPointSize (8);
      }

      // draw a point
      glColor3fv (cc->getAverage ());
      glBegin (GL_POINTS);
      glVertex2f (pt[0], pt[1]);
      glEnd ();
    }

    // do not draw the children if the node is marked as collapsed
    if ((*hPlot)[n->getCID ()].state &= TreePlotter::CLPS)
    {
      return;
    }

    // draw children
    if (l != NULL)
    {
      drawMSTNode (l, selected);
    }

    if (r != NULL)
    {
      drawMSTNode (r, selected);
    }
  }
public:
  SegmentationTree (const unsigned int textureID)
    : refTexture (textureID),
      hTree (NULL),
      hPlot (NULL),
      componentSet (NULL),
      selectedNode (NULL),
      selectedComponent (0),
      layoutStyle (RADIAL),
      drawStyle (BLOCK),
      split (0.01f)
  {
    createStipplePattern ();
  }

  ~SegmentationTree ()
  {
    if (hTree != NULL) {delete hTree; hTree = NULL;}
    if (hPlot != NULL) {delete hPlot; hPlot = NULL;}
    if (dcm != NULL) {delete dcm;   dcm = NULL;}
    if (componentSet != NULL) {delete componentSet;  componentSet = NULL;}
  }
  /**
   * set the point set to be used by the segmentation tree
   */
  void setMSTHierarchy (Tree *mst)
  {
    hTree = mst;
  }

  void setDataComponentMap (DataComponentMap *newdcm)
  {
    dcm = newdcm;
  }

  void setComponentSet (ComponentSet *newComponentSet)
  {
    componentSet = newComponentSet;
  }

  const ComponentSet* getComponentSet () const
  {
    return componentSet;
  }

  void buildHierarchyPlot (TreePlotter::PlotMethod method)
  {
    TreePlotter::TreePlotIterator  it;
    float  x[2], y[2];//, mn[2], mx[2];

    if (hTree == NULL)
      return;

  //  if (hPlot != NULL)  {delete hPlot; hPlot = NULL;}

    if (hPlot == NULL)
      hPlot = new TreePlotter::TreePlot ();

  //  TreePlotter::getTreeCoords_Basic (hTree, hPlot);
    TreePlotter::getTreeCoords (hTree, hPlot, method);

    // map the coords to the view
    tp_rng[0] = tp_rng[1] = FLT_MAX;
    tp_rng[2] = tp_rng[3] = -FLT_MAX;
    for (it = hPlot->begin (); it != hPlot->end (); it++)
    {
    //  x = ((TreePlotter::DParm)it->second).pt[0];
      x[0] = ((TreePlotter::DParm)it->second).xr[0];
      x[1] = ((TreePlotter::DParm)it->second).xr[1];
      y[0] = ((TreePlotter::DParm)it->second).pt[1];
      y[1] = ((TreePlotter::DParm)it->second).pt[1]+1.0f;

      if (!(it->second.state & TreePlotter::VALID))
        continue;

      if (it->second.state & TreePlotter::CLPS)
        continue;

      if (tp_rng[0]>x[0])  {tp_rng[0]=x[0];}
      if (tp_rng[1]>y[0])  {tp_rng[1]=y[0];}
      if (tp_rng[2]<x[1])  {tp_rng[2]=x[1];}
      if (tp_rng[3]<y[1])  {tp_rng[3]=y[1];}
    }
#if 1
    // unitize
    for (it = hPlot->begin (); it != hPlot->end (); it++)
    {
      if (!(it->second.state & TreePlotter::VALID))
        continue;

      x[0] = ((TreePlotter::DParm)it->second).pt[0];
      y[0] = ((TreePlotter::DParm)it->second).pt[1];

      x[0] = (x[0]-tp_rng[0])/(tp_rng[2]-tp_rng[0]);
      y[0] = (y[0]-tp_rng[1])/(tp_rng[3]-tp_rng[1]);

      (*hPlot)[it->first].pt[0] = x[0];
      (*hPlot)[it->first].pt[1] = y[0];//= TreePlotter::Point (x, y);

      x[0] = ((TreePlotter::DParm)it->second).xr[0];
      x[1] = ((TreePlotter::DParm)it->second).xr[1];

      x[0] = (x[0]-tp_rng[0])/(tp_rng[2]-tp_rng[0]);
      x[1] = (x[1]-tp_rng[0])/(tp_rng[2]-tp_rng[0]);

    //  y = 1.0f-y;
/*
      if (layoutStyle == RADIAL)
      {
        float  r, th;

        // make it radial
        r = y;
        // rotate the point
        th = (float)(x*2.0*M_PI); // 2.0 * ;
        y = (r) * (cos (th) - sin (th));
        x = (r) * (cos (th) + sin (th));

        th = (float)(xr[0]*2.0*M_PI);  xr[0] = (r) * (cos (th) + sin (th));
        th = (float)(xr[1]*2.0*M_PI);  xr[1] = (r) * (cos (th) + sin (th));
      }
*/
      // copy back
      (*hPlot)[it->first].xr[0] = x[0];
      (*hPlot)[it->first].xr[1] = x[1];
    }
#endif
    tp_rng[0] = tp_rng[1] = FLT_MAX;
    tp_rng[2] = tp_rng[3] = -FLT_MAX;
    for (it = hPlot->begin (); it != hPlot->end (); it++)
    {
    //  x = ((TreePlotter::DParm)it->second).pt[0];
      x[0] = ((TreePlotter::DParm)it->second).xr[0];
      x[1] = ((TreePlotter::DParm)it->second).xr[1];
      y[0] = ((TreePlotter::DParm)it->second).pt[1];
    //  y[1] = ((TreePlotter::DParm)it->second).pt[1];

      if (!(it->second.state & TreePlotter::VALID))
        continue;

      if (it->second.state & TreePlotter::CLPS)
        continue;

      if (tp_rng[0]>x[0])  {tp_rng[0]=x[0];}
      if (tp_rng[1]>y[0])  {tp_rng[1]=y[0];}
      if (tp_rng[2]<x[1])  {tp_rng[2]=x[1];}
      if (tp_rng[3]<y[0])  {tp_rng[3]=y[0];}
    }
  }
  // drawing
  /**
   * setup the 3d projection
   * @param camera camera object specifying the view controls
   */
  void setupView (const GLUtils::Camera& camera)
  {
  //  float  vp[4];

  //  glGetFloatv (GL_VIEWPORT, vp);

    // get the viewport size
    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    glLoadIdentity ();
  //  gluPerspective (45.0, vp[2]/vp[3], 1, 100);
    glOrtho (-2, 2, -2, 2, 0, 1);

    glMatrixMode (GL_MODELVIEW);
    glPushMatrix ();
    glLoadIdentity ();
  //  camera.rotate ();
  }
  /**
   * pop the modelview and projection matrices pushed onto the stack in setup3DView
   */
  void cleanupView ()
  {
    glMatrixMode (GL_PROJECTION);
    glPopMatrix ();

    glMatrixMode (GL_MODELVIEW);
    glPopMatrix ();
  }
  /**
   * draw the reference texture
   */
  void drawTexture ()
  {
    float  vp[4];
    int    w, h;

    if (glIsTexture (refTexture) == GL_TRUE)
    {
      glGetFloatv (GL_VIEWPORT, vp);

      glMatrixMode (GL_PROJECTION);  glPushMatrix ();  glLoadIdentity ();  glOrtho (0, vp[2], vp[3], 0, 0, 1);
      glMatrixMode (GL_MODELVIEW);  glPushMatrix ();  glLoadIdentity ();

      glEnable (GL_TEXTURE_2D);
      glBindTexture (GL_TEXTURE_2D, refTexture);
      glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

      glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
      glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

      glBegin (GL_QUADS);
      glTexCoord2i (0, 1);  glVertex2i (0, h);
      glTexCoord2i (0, 0);  glVertex2i (0, 0);
      glTexCoord2i (1, 0);  glVertex2i (w, 0);
      glTexCoord2i (1, 1);  glVertex2i (w, h);
      glEnd ();

      glDisable (GL_TEXTURE_2D);

      glMatrixMode (GL_PROJECTION);  glPopMatrix ();
      glMatrixMode (GL_MODELVIEW);  glPopMatrix ();
    }
  }

  /**
   * draw the MST Hierarchy
   */
  void drawMSTTree ()
  {
    TreePlotter::TreePlotIterator  it;
    TreePlotter::Point        mn, mx, fr;

    if ((hPlot == NULL) || (hPlot->empty ()))
      return;

    {
    // background
    glMatrixMode (GL_PROJECTION);  glPushMatrix ();  glLoadIdentity ();  glOrtho (0, 1, 1, 0, 0, 1);
    glMatrixMode (GL_MODELVIEW);  glPushMatrix ();  glLoadIdentity ();

    glBegin (GL_QUADS);
    glColor3f (1,1,1);  glVertex2f (0,0);  glVertex2f (1, 0);
    glColor3f (1,1,1);  glVertex2f (1,1);  glVertex2f (0, 1);
    glEnd ();

    glMatrixMode (GL_PROJECTION);  glPopMatrix ();
    glMatrixMode (GL_MODELVIEW);  glPopMatrix ();

    // find the min and max of the tree
    mn[0] = tp_rng[0];    mn[1] = tp_rng[1];
    mx[0] = tp_rng[2];    mx[1] = tp_rng[3];

    if (layoutStyle == RADIAL)
    {
      mn[0]=mn[1]=-M_PI_2;
      mx[0]=mx[1]= M_PI_2;
    }
    fr[0] = mx[0]-mn[0];  fr[1] = mx[1]-mn[1];
    fr[0] = 0.1f;      fr[1] = 0.1f;

    // setup the projection
    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    glLoadIdentity ();
    glOrtho (mn[0]-fr[0], mx[0]+fr[0], mx[1]+fr[1], mn[1]-fr[1], 0, 1);

    glMatrixMode (GL_MODELVIEW);
    glPushMatrix ();
    glLoadIdentity ();//  glTranslatef ((zm/2.0f), (zm/2.0f), 0.0f);

  //  camera->translate ();
    }
    glLineWidth (2);
    glPointSize (8);
    glEnable (GL_POINT_SMOOTH);
  //  glEnable (GL_ALPHA_TEST);
  //  glAlphaFunc (GL_GREATER, 0.5f);

    if (layoutStyle == RADIAL)
    {
      drawCircleTree ();
    }

    // drwa the nodes
    glEnable (GL_POLYGON_STIPPLE);
    drawMSTNode (hTree->root (), selectedNode == hTree->root ());
    glDisable (GL_POLYGON_STIPPLE);
/*
    else
    {  // old method
    //  for (i=0; i<hTree->length (); i++)
      for (it = hPlot->begin (); it != hPlot->end (); it++)
      {
        if (!(it->second.state & TreePlotter::VALID))
          continue;

        pt1[0] = it->second.pt[0];
        pt1[1] = it->second.pt[1];

        cc = componentSet->find (it->first);
        n = hTree->find (it->first);
        // draw an edge to the parent
  #if 1
        if (n->parent () != NULL)
        {
          cp = componentSet->find (n->parent ()->getCID ());

          // find p in the tree
          pt2[0] = (*hPlot)[n->parent ()->getCID ()].pt[0];
          pt2[1] = (*hPlot)[n->parent ()->getCID ()].pt[1];

          // draw the connection
          glBegin (GL_LINES);
          glColor3fv (cc->getAverage ());  glVertex2fv (pt1);
          glColor3fv (cp->getAverage ());  glVertex2fv (pt2);
          glEnd ();
        }
  #endif
        // draw a triangle between the parent and left and right children
  #if 0
        if ((n->left () != NULL) && (n->right () != NULL))
        {
          cl = componentSet->find (n->left ()->getCID ());
          cr = componentSet->find (n->right ()->getCID ());

          glBegin (GL_TRIANGLES);
          glColor3fv (cc->getAverage ());  glVertex2fv (pt1);
          glColor3fv (cl->getAverage ());  glVertex2fv ((*hPlot)[n->left ()->getCID ()].pt);
          glColor3fv (cr->getAverage ());  glVertex2fv ((*hPlot)[n->right ()->getCID ()].pt);
          glEnd ();
        }
  #endif

        // draw the node
  #if 0
        if (cc->hasTexture ())
        {
          pt[0] = pt1[0];
          pt[1] = pt1[1];
          s = 0.1f * (1.0f-pow (pt[1], 4.0f));

          // draw a label
          glEnable (GL_TEXTURE_2D);
          glBindTexture (GL_TEXTURE_2D, cc->getTexID ());//ccTex[ccMap[c]]);
          glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

          glBegin (GL_QUADS);
          glTexCoord2f (0,0);  glVertex2f (pt[0]-s, pt[1]-s);
          glTexCoord2f (0,1);  glVertex2f (pt[0]-s, pt[1]+s);
          glTexCoord2f (1,1);  glVertex2f (pt[0]+s, pt[1]+s);
          glTexCoord2f (1,0);  glVertex2f (pt[0]+s, pt[1]-s);
          glEnd ();
          glDisable (GL_TEXTURE_2D);
        }
        else
  #endif
        {
          if (isSelected (it->first))
          {
            glPointSize (12);
          //  glColor3f (1,1,1);
            glColor3f (0,0,0);
            glBegin (GL_POINTS);
            glVertex2fv (pt1);
            glEnd ();
            glPointSize (8);
          }

          glColor3fv (cc->getAverage ());

          // draw a point
          glBegin (GL_POINTS);
          glVertex2fv (pt1);// (pt1[0], pt1[1]);
          glEnd ();
        }
      }
    }

*/
  //  glDisable (GL_ALPHA_TEST);
    glDisable (GL_POINT_SMOOTH);
    glPointSize (1);
    glLineWidth (1);

    // cleanup
    glMatrixMode (GL_PROJECTION);  glPopMatrix ();
    glMatrixMode (GL_MODELVIEW);  glPopMatrix ();
  }

  /**
   * draw the set of leaves as an image
   */
  void drawLeafComponents ()
  {
    ComponentSet::Component *c, *cp;
    const Tree::Node *n, *p;
    unsigned int i;
    const float *co;
    float o[2];  // offset
    float col[4], sel_col[4];

  //  if (!hasComponentImages ())
  //    return;

    col[0] = col[1] = col[2] = 1.0f; col[3] = 1.0f;
    sel_col[0] = sel_col[1] = sel_col[2] = sel_col[3] = 1.0f;

    glEnable (GL_DEPTH_TEST);
    glEnable (GL_ALPHA_TEST);
    glAlphaFunc (GL_GREATER, 0.5f);

    glEnable (GL_TEXTURE_2D);

    for (i=0; i<hTree->length (); i++)
    {
      if ((n = hTree->get (i)) == NULL)
        continue;

      if (!n->isLeaf ())
        continue;

      if ((c = componentSet->find (n->getCID ())) == NULL)
        continue;

      if (c->hasTexture () == false)
        continue;

      // get the offset from the parent down
      for (o[0]=0.0f, o[1]=0.0f, p=n; p!=NULL; p=p->parent ())
      {
        // if not node is selected, do not separate the children
        // by removing their contribution to the vector
        if (isSelected (p->getCID ()))
        {
          o[0] = 0.0f;
          o[1] = 0.0f;
        }

        // if the component is collapsed, do not separate the children
        if ((*hPlot)[p->getCID ()].state &= TreePlotter::CLPS)
        {
          o[0] = 0.0f;
          o[1] = 0.0f;
        }

        if ((cp = componentSet->find (p->getCID ())) == NULL)
          continue;

        co = cp->getOffset ();

        o[0] += co[0]*split;//ccOffs[(ccMap[p->getCID ()]*2)+0]*split;
        o[1] += co[1]*split;//ccOffs[(ccMap[p->getCID ()]*2)+1]*split;
      }


      glPushMatrix ();
      glTranslatef (o[0], o[1], 0.0f);//, o[1]);

      // draw a label
      glBindTexture (GL_TEXTURE_2D, c->getTexID ());//getComponentTexID (n->getCID ()));//ccTex[ccMap[c]]);
      glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    //  // setup the texture mode
    //  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);
    //  // get the alpha from the texture
    //  glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_REPLACE);
    //  glTexEnvi(GL_TEXTURE_ENV, GL_SRC0_ALPHA, GL_TEXTURE);
    //  glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_ALPHA, GL_SRC_ALPHA);
    //  // use the colour from the texture unless we are a selected node,
    //  // then use the node colour
    //  glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB, GL_MODULATE);
    ////  glTexEnvi (GL_TEXTURE_ENV, GL_SRC0_RGB, GL_CONSTANT);
    ////  glTexEnvi (GL_TEXTURE_ENV, GL_OPERAND0_RGB, GL_SRC_COLOR);

    //  if (isSelected (n->getCID ()))
    //  {
    //    glTexEnvi (GL_TEXTURE_ENV, GL_SRC0_RGB, GL_TEXTURE);
    //    glTexEnvi (GL_TEXTURE_ENV, GL_OPERAND0_RGB, GL_CONSTANT);

    //  //  glColor3f (1,0,0);
    //    cp = componentSet->find (getSelectedComponent ());
    //    memcpy (sel_col, cp->getAverage (), 3*sizeof (float));
    //    glColor3fv (sel_col);
    //    glTexEnvfv (GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, sel_col);
    //  }
    //  else
    //  {
    //    glTexEnvi (GL_TEXTURE_ENV, GL_SRC0_RGB, GL_TEXTURE);
    //    glTexEnvi (GL_TEXTURE_ENV, GL_OPERAND0_RGB, GL_SRC_COLOR);
    //  //  glColor3f (0,0,0);
    //    glColor4fv (col);
    //    glTexEnvfv (GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, col);
    //  }

      glBegin (GL_QUADS);
      glTexCoord2f (0,0);  glVertex2f (-1.0f, -1.0f);
      glTexCoord2f (0,1);  glVertex2f (-1.0f,  1.0f);
      glTexCoord2f (1,1);  glVertex2f ( 1.0f,  1.0f);
      glTexCoord2f (1,0);  glVertex2f ( 1.0f, -1.0f);
      glEnd ();

      glPopMatrix ();
    }

    glDisable (GL_TEXTURE_2D);
    glDisable (GL_ALPHA_TEST);
    glDisable (GL_DEPTH_TEST);
  }

  // draw points
  void drawPoints ()
  {
    int    x, y;
    int    w, h;
    float  px, py, tx, ty;

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable (GL_TEXTURE_2D);
    glBindTexture (GL_TEXTURE_2D, refTexture);
    glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    glBegin (GL_QUADS);

    for (y=0; y<h; y++)
    {
      for (x=0; x<w; x++)
      {
        if (!isSelected (dcm->getDataComponent ((y*w)+x)))
          glColor4f (1,1,1,0.1);
        else
          glColor4f (1,1,1,1);

        tx=((float)x/(float)w);  ty=((float)y/(float)h);
        glTexCoord2f (tx, ty);  glVertex2f ((tx*2.0f)-1.0f, (ty*2.0f)-1.0f);
        tx=((float)(x+1)/(float)w);  ty=((float)y/(float)h);
        glTexCoord2f (tx, ty);  glVertex2f ((tx*2.0f)-1.0f, (ty*2.0f)-1.0f);
        tx=((float)(x+1)/(float)w);  ty=((float)(y+1)/(float)h);
        glTexCoord2f (tx, ty);  glVertex2f ((tx*2.0f)-1.0f, (ty*2.0f)-1.0f);
        tx=((float)x/(float)w);  ty=((float)(y+1)/(float)h);
        glTexCoord2f (tx, ty);  glVertex2f ((tx*2.0f)-1.0f, (ty*2.0f)-1.0f);
      }
    }
    glEnd ();

    glDisable (GL_TEXTURE_2D);
    glDisable (GL_BLEND);
  }

  Tree *getMSTHierarchy ()
  {
    return hTree;
  }

  DataComponentMap *getDataComponentMap () const
  {
    return dcm;
  }

  const unsigned int getMSTHierarchyCount () const
  {
  //  return mstHierarchyCnt;
    if (hTree != NULL)
      return hTree->length ();
    else
      return 0;
  }

  /**
   * indicates whether the node with the component id (or any of its ancestors)
   * is selected.
   * @return true if the component, or any ancestors, is selected.
   */
  const bool isSelected (const unsigned int CID) const
  {
    Tree::Node  *tc;

    if (selectedNode == NULL)
      return false;

  //  if (selectedNode->getCID () == CID)
  //    return true;
  //  else
  //    return false;

    tc = hTree->find (CID);

    if (tc == NULL)
      return false;

    if (!tc->isAncestor (selectedNode))
      return false;

    return true;
  }

  // interaction
  const unsigned int getSelectedComponent () const
  {
    return selectedComponent;
  }

  void setSelectedComponent (const unsigned int newSelectedComponent)
  {
    selectedComponent = newSelectedComponent;

    selectedNode = hTree->find (selectedComponent);
  }

  void selectParentComponent (const unsigned int childComponent)
  {
  //  if (mstHierarchy == NULL)
  //    return;

  //  selectedComponent = mstHierarchy[childComponent];
    if (hTree  == NULL)
      return;

    Tree::Node  *n;

    n = hTree->find (childComponent);

    if ((n == NULL) || (n->parent () == NULL))
      return;

    setSelectedComponent (n->parent ()->getCID ());
  }

  void selectLeftComponent (const unsigned int parentComponent)
  {
    if (hTree == NULL)
      return;

    Tree::Node  *n;

    n = hTree->find (parentComponent);

    if ((n == NULL) || (n->left () == NULL))
      return;

    setSelectedComponent (n->left ()->getCID ());
  }

  void selectRightComponent (const unsigned int parentComponent)
  {
    if (hTree == NULL)
      return;

    Tree::Node  *n;

    n = hTree->find (parentComponent);

    if ((n == NULL) || (n->right () == NULL))
      return;

    setSelectedComponent (n->right ()->getCID ());
  }

  void setLayoutStyle (LayoutStyle newLayoutStyle)
  {
    layoutStyle = newLayoutStyle;
  }

  const LayoutStyle getLayoutStyle () const
  {
    return layoutStyle;
  }

  void setDrawStyle (DrawStyle newDrawStyle)
  {
    drawStyle = newDrawStyle;
  }

  const DrawStyle getDrawStyle () const
  {
    return drawStyle;
  }

  const float getSeparation () const
  {
    return split;
  }

  void setSeparation (const float newSeparation)
  {
    split = newSeparation;
  }

  void collapseSelected ()
  {
    if (selectedNode == NULL)
      return;

    if ((*hPlot)[selectedNode->getCID ()].state &= TreePlotter::CLPS)
      (*hPlot)[selectedNode->getCID ()].state ^= TreePlotter::CLPS;
    else
      (*hPlot)[selectedNode->getCID ()].state |= TreePlotter::CLPS;
  }
};

