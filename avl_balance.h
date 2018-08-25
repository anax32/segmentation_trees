class AVLBalance
{
protected:
  Tree *tree;
  std::stack<Tree::Node*> cn;
//  std::map<unsigned int, bool>  cnp;

  int factor (Tree::Node *n)
  {
    int  lh, rh;

    if (n->left () != NULL)
      lh = n->left ()->height ();
    else
      lh = 0;

    if (n->right () != NULL)
      rh = n->right ()->height ();
    else
      rh = 0;

    return (lh-rh);
  }

  void rotateLeft (Tree::Node *n)
  {
    Tree::Node *r = NULL;
    Tree::Node *rl = NULL;
    Tree::Node *np = NULL;

    if (n->right () == NULL)
      return;
/*
    a
     \                   b
      b        ->       / \
       \         a   c
        c
*/
    np = n->parent ();
    r = n->right ();
    rl = n->right ()->left ();

    n->setRight (rl);
    r->setLeft (n);
    r->setParent (np);

    // change the pointer in the parent
    if (np != NULL)
    {
      if (np->left () == n)
        np->setLeft (r);
      else if (np->right () == n)
        np->setRight (r);
      else
        n=n;  // error!
    }
  }
  void rotateRight (Tree::Node *n)
  {
    Tree::Node  *l = NULL;
    Tree::Node  *lr = NULL;
    Tree::Node  *np = NULL;

    if (n->left () == NULL)
      return;
/*
      c
       /        b
        b      ->     / \
       /          a   c
      a
*/
    np = n->parent ();
    l = n->left ();
    lr = n->left ()->right ();

    // change the pointer in the parent
    n->setLeft (lr);
    l->setParent (np);
    l->setRight (n);

    if (np != NULL)
    {
      if (np->left () == n)
        np->setLeft (l);
      else if (np->right () == n)
        np->setRight (l);
      else
        n=n;  // error!
    }
  }

  void balanceNode (Tree::Node *n)
  {
    int  nf, rf, lf;

  //  fprintf (stdout, "CID: %04i/%04i\n\0", n->getCID (), tree->length ());

    // balance the children
    if ((n->right () != NULL))// && (cnp[n->right ()->getCID ()] == false))
    {
    //  balanceNode (n->right ());

      rf = factor (n->right ());

      if (abs (rf) > 1)
      {
        cn.push (n->right ());
      //  cnp[n->right ()->getCID ()] = true;
        return;
      }
    }

    if ((n->left () != NULL))// && (cnp[n->left ()->getCID ()] == false))
    {
      //  balanceNode (n->left ());
      lf = factor (n->left ());

      if (abs (lf) > 1)
      {
        cn.push (n->left ());
      //  cnp[n->left ()->getCID ()] = true;
        return;
      }
    }

    n = cn.top ();
  //  cnp[n->getCID ()] = true;

    // balance this node
    nf = factor (n);

    if (abs (nf) < 2)
    {
      cn.pop ();
      return;
    }

    // need balancing
    if (nf < 0)    // tree is right heavy
    {
      // check if the right subtree is left heavy
      rf = factor (n->right ());

      if (rf > 0)    // tree's right subtree is left heavy
      {
      // Perform Double Left rotation (right-left rotation)
        rotateRight (n->right ());
        rotateLeft (n);
      }
      else
      {
      // Perform Single Left rotation
        rotateLeft (n);
      }
    }
    else if (nf > 0) // tree is left heavy
    {
      lf = factor (n->left ());

      if (lf < 0)  //  tree's left subtree is right heavy
      {
      //  Perform Double Right rotation (left-right rotation)
        rotateLeft (n->left ());
        rotateRight (n);
      //  rotateRight (n);
      }
      else
      {
      //  Perform Single Right rotation
        rotateRight (n);
      }
    }

    n = n->parent ();

    // reset
    cn.pop ();

    if (n != NULL)
    {
      cn.push (n);
    //  cnp[n->getCID ()] = false;
    }
  //  if (abs (factor (n)) < 2)
  //  {
  //    cn.push (n);
  //    cnp[n->getCID ()] = false;
  //  }
  }

public:
  AVLBalance (Tree *unbalanced)
    : tree (unbalanced)
  {}

  ~AVLBalance ()
  {}

  void balance ()
  {
    unsigned int  t;

    fprintf (stdout, "Balancing tree...\0");
    t = clock ();
    balanceNode (const_cast<Tree::Node*>(tree->root ()));

    // find the root
  //  for (i=0; i<tree->length (); i++)
  //  {
  //    if (tree->getNodeIdx (i)->parent () == NULL)
  //      tree->setRootID (i);
  //  }

    fprintf (stdout, "\rBalanced tree in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
  }

  bool next ()
  {
    if (cn.empty ())
    {
      cn.push (const_cast<Tree::Node*>(tree->root ()));
    //  cnp[tree->root ()->getCID ()] = true;
    }

    balanceNode (cn.top ());

  //  if (cn.top ()->parent () == NULL)
    if (cn.empty ())
    {
      tree->setRoot ();
      return false;
    }

    return true;
  }

  Tree::Node *getCN ()
  {
    return cn.top ();
  }
};
