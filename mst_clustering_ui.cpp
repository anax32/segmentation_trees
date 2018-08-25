#pragma comment(lib, "../../Image/PNG/libpng13.lib")  // libpng

#include <iostream>

#define NOMINMAX
#define OEMRESOURCE
#include <windows.h>
#include <direct.h>
#include <assert.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "../../Utils/Patterns/Persistable.h"
#include "../../Utils/Structures/ParameterSet.h"
#include "../../Utils/Callbacks/Callback.h"

#include "../../Utils/Win/Utils.h"  // console
#include "../../Utils/Win/Messages.h"
#include "../../Utils/Win/Console.h"
#include "../../Image/pngreader.h"  // libpng
#include "../../Utils/ImageUtils.h"
//#include "../../DataStructures/Queue.h"
#include "../../Utils/GraphUtils.h"
#include "../../Utils/ColourUtils.h"

// gl
#define _GLUTILS_DEF_LIBS_
#define _GLUTILS_USE_GLEW_
#include "../../Utils/GL/GLUtils.h"
#include "../../Utils/GL/Camera.h"

#include "../DataComponentMap.h"
#include "../Tree.h"
#include "../Components.h"
#include "../SegmentationTree.h"
#include "../MSTree.h"
#include "../AVLBalance.h"

//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\cube.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\baboon_128.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\baboon_256.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\baboon_512.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\glacier1_512.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\glacier1_256.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\glacier1_128.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\glacier1_64.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\bearglacier_sm.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\lena.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\lena_512.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\cameraman_512.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\kutztown_1024.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\kutztown_512.png\0"
//#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\imgs\\kutztown_256.png\0"

#define CACHE_DIR  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\cache\0"

#define PRUNE_DEPTH  1024
#undef NO_PRUNE
/*===============================================
 Construct a graph from a gl texture
===============================================*/
void constructPointSet (unsigned int texture, float **pointSet, unsigned int& pointCount, unsigned int& pointDims,
            const bool uniques, DataComponentMap *dcm)
{
  // build a point set from a thresholded texture
  unsigned int  x, y, i, t, pi;
  unsigned int  w, h, s;
  float      *img;
  float      *p, px[5];

  fprintf (stdout, "Constructing point set...\0");
  t = clock ();

  // first construct a point set of the texture
  glBindTexture (GL_TEXTURE_2D, texture);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, (GLint*)&w);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, (GLint*)&h);
  s=w*h;

  // allocate space
  img = new float[s*3];
  glGetTexImage (GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, img);

  // allocate space
  p = new float[s*5];
  memset (p, 0, s*5*sizeof (float));

  pi = 0;

  for (y=0;y<h;y++)
  {
    for (x=0;x<w;x++)
    {
      dcm->setDataComponent ((y*w)+x, (y*w)+x);

      px[0] = ((float)x/(float)w)*1.0f;  //0.0f;
      px[1] = ((float)y/(float)h)*1.0f;  //0.0f;
      px[2] = img[(y*w*3)+(x*3)+0];
      px[3] = img[(y*w*3)+(x*3)+1];
      px[4] = img[(y*w*3)+(x*3)+2];

      if (uniques)
      {
        // check if this point already exists in the data
        for (i=0; i<pi; i++)
        {
          if ((p[(i*5)+0]==px[0]) && (p[(i*5)+1]==px[1]) && (p[(i*5)+2]==px[2]) &&
            (p[(i*5)+3]==px[3]) && (p[(i*5)+4]==px[4]))
            break;
        }

        // yes
        if (i<pi)
        {
          dcm->setDataComponent ((y*w)+x, i);
          continue;
        }
        // no, continue...
      }

      // add the item to the pointlist
      i = pi*5;
      p[i+0] = px[0];
      p[i+1] = px[1];
      p[i+2] = px[2];
      p[i+3] = px[3];
      p[i+4] = px[4];
      pi++;
    }
  }

  delete[] img;
  img = NULL;

  *pointSet = p;
  pointCount = pi;
  pointDims = 5;

  fprintf (stdout, "\rConstructed point set in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
}

float vertexWeightd (const float *v1, const  float *v2)
{
  const float    pdwt = 8.0f;
  float  rgb_p1[3], rgb_p2[3];
  float  yuv_p1[3], yuv_p2[3];
  float  pd, lc, cd;
//  float  cc[3];

  rgb_p1[0] = v1[2];
  rgb_p1[1] = v1[3];
  rgb_p1[2] = v1[4];

  rgb_p2[0] = v2[2];
  rgb_p2[1] = v2[3];
  rgb_p2[2] = v2[4];
#if 0
  ColourUtils::rgb2yuv (rgb_p1, yuv_p1);
  ColourUtils::rgb2yuv (rgb_p2, yuv_p2);

  // yuv distance
  lc = (yuv_p1[0]-yuv_p2[0])*(yuv_p1[0]-yuv_p2[0]);
  lc += (yuv_p1[1]-yuv_p2[1])*(yuv_p1[1]-yuv_p2[1]);
  lc += (yuv_p1[2]-yuv_p2[2])*(yuv_p1[2]-yuv_p2[2]);

  if (lc > 0.0f)
    lc = sqrt (lc);
#else
  lc = 0;
#endif

  // colour contrast
//  cc[0] = abs (rgb_p1[0]/(rgb_p1[0]-rgb_p1[1])) - abs (rgb_p2[0]/(rgb_p2[0]-rgb_p2[1]));
//  cc[1] = abs (rgb_p1[1]/(rgb_p1[1]-rgb_p1[0])) - abs (rgb_p2[1]/(rgb_p2[1]-rgb_p2[0]));

  // physical distance
  pd =  (v1[0]*pdwt)-(v2[0]*pdwt);
  pd += (v1[1]*pdwt)-(v2[1]*pdwt);

  if (pd > 0.0f)
    pd = sqrt (pd);

  // rgb distance
  cd = (rgb_p1[0]-rgb_p2[0])*(rgb_p1[0]-rgb_p2[0]);
  cd += (rgb_p1[1]-rgb_p2[1])*(rgb_p1[1]-rgb_p2[1]);
  cd += (rgb_p1[2]-rgb_p2[2])*(rgb_p1[2]-rgb_p2[2]);

  if (cd > 0.0f)
    cd = sqrt (cd);

  // get the distance
//  return lc; //sqrt (/*x*x+y*y+*/lc);//+cc[0]*cc[0]+cc[1]*cc[1]);
  return pd + cd + lc;
}

float vertexWeight (const float *vertices, const unsigned int v1, const unsigned int v2)
{
#if 0
  float  x,y;
  float  rgb_p1[3], rgb_p2[3];
  float  yuv_p1[3], yuv_p2[3];
  float  lc, cd;
//  float  cc[3];

  x = (vertices[(v1*5)+0]-vertices[(v2*5)+0])*2.0f;
  y = (vertices[(v1*5)+1]-vertices[(v2*5)+1])*2.0f;
  rgb_p1[0] = vertices[(v1*5)+2];
  rgb_p1[1] = vertices[(v1*5)+3];
  rgb_p1[2] = vertices[(v1*5)+4];

  rgb_p2[0] = vertices[(v2*5)+2];
  rgb_p2[1] = vertices[(v2*5)+3];
  rgb_p2[2] = vertices[(v2*5)+4];

  ColourUtils::rgb2yuv (rgb_p1, yuv_p1);
  ColourUtils::rgb2yuv (rgb_p2, yuv_p2);

  // yuv distance
  lc = (yuv_p1[0]-yuv_p2[0])*(yuv_p1[0]-yuv_p2[0]);
  lc += (yuv_p1[1]-yuv_p2[1])*(yuv_p1[1]-yuv_p2[1]);
  lc += (yuv_p1[2]-yuv_p2[2])*(yuv_p1[2]-yuv_p2[2]);

  if (lc > 0.0f)
    lc = sqrt (lc);

  // colour contrast
//  cc[0] = abs (rgb_p1[0]/(rgb_p1[0]-rgb_p1[1])) - abs (rgb_p2[0]/(rgb_p2[0]-rgb_p2[1]));
//  cc[1] = abs (rgb_p1[1]/(rgb_p1[1]-rgb_p1[0])) - abs (rgb_p2[1]/(rgb_p2[1]-rgb_p2[0]));

  // rgb distance
  cd = (rgb_p1[0]-rgb_p2[0])*(rgb_p1[0]-rgb_p2[0]);
  cd += (rgb_p1[1]-rgb_p2[1])*(rgb_p1[1]-rgb_p2[1]);
  cd += (rgb_p1[2]-rgb_p2[2])*(rgb_p1[2]-rgb_p2[2]);

  if (cd > 0.0f)
    cd = sqrt (cd);

  // get the distance
//  return lc; //sqrt (/*x*x+y*y+*/lc);//+cc[0]*cc[0]+cc[1]*cc[1]);
  return cd;// + lc;
#else
  return vertexWeightd (&vertices[v1*5], &vertices[v2*5]);
#endif
}

/*===============================================
 Component info
===============================================*/
ComponentSet* buildComponents (const char *cacheDir, const unsigned int refTexture, Tree *tree, DataComponentMap *dcm)
{
  ComponentSet::Component  *cc;
  ComponentSet      *cs;
//  unsigned char  *msk, *img;
  unsigned int  i, j, t;
  Tree::Node    *n;

  fprintf (stdout, "Building components...\0");
  t = clock ();

  // allocate the component object
  cs = new ComponentSet (tree->length (), 2);

  //// alloc a buffer for the alpha mask
  //glBindTexture (GL_TEXTURE_2D, refTexture);
  //glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
  //glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

  //msk = new unsigned char[w*h*4];
  //img = new unsigned char[w*h*3];

  //// get the reference image
  //glGetTexImage (GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_BYTE, img);

  fprintf (stdout, "\rLinking components to nodes...\0");

  // build the component list
  for (i=0, j=0; i<tree->length (); i++)
  {
    // get the node
    n = tree->get (i);

  //  if (n->depth () > minDepth)
  //  if (n->countLeaves () < minVol)
  //  if (!(n->isLeaf ()) /*&& (hTree->isParent (n))*/)
  //    continue;

    if (!n->isAncestor (tree->root ()))
      continue;

  //  if (!n->isLeaf ())
  //    continue;

    // only set the cids here to link the components to nodes
    cc = cs->get (j++);
    cc->setCID (n->getCID ());

#if 0
    char  buf[MAX_PATH];

    sprintf_s (buf, MAX_PATH, "%sccsimg\\%i-cc.png\0", cacheDir, n->getCID ());
    GLUtils::writeTextureToFile (buf, GL_RGBA);
#endif
  }

  fprintf (stdout, "\rLinked components to nodes in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);

  // cleanup
//  delete[] msk;
//  delete[] img;

  // next
  return cs;
}

void getNodeTabs (ComponentSet *cs, Tree *tree, const Tree::Node *n, const unsigned int tw, const unsigned int th)
{
//  const float  *or;

  // project the component into the unit square to get a square representation
//  or = cc->getOrigin ();
  ComponentSet::Component  *cc;
  unsigned int  i, ii, x, y;
  float      *tab, *tl, *tr;
  const float    *avg;

  cc = cs->find (n->getCID ());
  tab = new float[tw*th*4];
  memset (tab, 0, tw*th*4*sizeof (float));

  // build tab images for the components to serve as a small summary
  if (n->isLeaf ())
  {
    avg = cc->getAverage ();

    for (i=0, ii=0; i<tw*th; i++, ii+=4)
    {
      tab[ii+0] = avg[0];
      tab[ii+1] = avg[1];
      tab[ii+2] = avg[2];
      tab[ii+3] = 1.0f;
    }
  }
  else
  {
    tl = tr = NULL;

    // recurse
    if (n->left () != NULL)
    {
      getNodeTabs (cs, tree, n->left (), tw, th);

      tl = new float[tw*th*4];
      glBindTexture (GL_TEXTURE_2D, cs->find (n->left ()->getCID ())->getTabID ());
      glGetTexImage (GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, tl);
    }
    if (n->right () != NULL)
    {
      getNodeTabs (cs, tree, n->right (), tw, th);

      tr = new float[tw*th*4];
      glBindTexture (GL_TEXTURE_2D, cs->find (n->right ()->getCID ())->getTabID ());
      glGetTexImage (GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, tr);
    }

    // merge the two images
    if ((tr == NULL) && (tl == NULL))
    {
      return;
    }
    else if ((tr != NULL) && (tl == NULL))
    {
      memcpy (tab, tr, tw*th*4*sizeof (float));
    }
    else if ((tl != NULL) && (tr == NULL))
    {
      memcpy (tab, tl, tw*th*4*sizeof (float));
    }
    else if ((tl != NULL) && (tr != NULL))
    {
      // FIXME: mix the two images by rendering both as quads to an
      //      fbo (dst=tabid) and use alpha blending to mix the images
      float      *src, a;
      unsigned int  xls, xrs, xmid;
      float      lrgba[4], rrgba[4];

      float      ns, ls, rs;

      ns = n->size ();
      ls = n->left ()->size ();
    //  rs = n->right ()->size ();

      xmid = ((float)ls/(float)ns)*tw;

      for (y=0; y<th; y++)
      {
        for (x=0; x<tw; x++)
        {
#if 1
          if (x < xmid)
          {
            src = tl;
            xls = ((float)x/(float)xmid)*tw;
          }
          else
          {
            src = tr;
            xls = ((float)(x-xmid)/(float)(tw-xmid))*tw;
          }

          tab[(y*tw*4)+(x*4)+0] = src[(y*tw*4)+(xls*4)+0];
          tab[(y*tw*4)+(x*4)+1] = src[(y*tw*4)+(xls*4)+1];
          tab[(y*tw*4)+(x*4)+2] = src[(y*tw*4)+(xls*4)+2];
          tab[(y*tw*4)+(x*4)+3] = src[(y*tw*4)+(xls*4)+3];
#else
        //  a = (float)x/(float)(tw-1);
          if (x < xmid)
            a = ((float)(x)/(float)(xmid))*((float)xmid/(float)tw);
          else
            a = (float)(x)/(float)(tw);

          lrgba[0] = tl[(y*tw*4)+(x*4)+0];
          lrgba[1] = tl[(y*tw*4)+(x*4)+1];
          lrgba[2] = tl[(y*tw*4)+(x*4)+2];
          lrgba[3] = tl[(y*tw*4)+(x*4)+3];

          rrgba[0] = tr[(y*tw*4)+(x*4)+0];
          rrgba[1] = tr[(y*tw*4)+(x*4)+1];
          rrgba[2] = tr[(y*tw*4)+(x*4)+2];
          rrgba[3] = tr[(y*tw*4)+(x*4)+3];

          tab[(y*tw*4)+(x*4)+0] = lrgba[0] + (a * (rrgba[0]-lrgba[0]));
          tab[(y*tw*4)+(x*4)+1] = lrgba[1] + (a * (rrgba[1]-lrgba[1]));
          tab[(y*tw*4)+(x*4)+2] = lrgba[2] + (a * (rrgba[2]-lrgba[2]));
          tab[(y*tw*4)+(x*4)+3] = lrgba[3] + (a * (rrgba[3]-lrgba[3]));
#endif
        }
      }
    }
    else
    {
      return;
    }
  }

  // upload
  cc->setTabData (tw, th, tab, GL_RGBA, GL_FLOAT);

#if 0
  glBindTexture (GL_TEXTURE_2D, cc->getTabID ());
  GLUtils::writeTextureToFile ("test.png\0", GL_RGBA);
#endif
  // cleanup
  delete[] tab;
}

void getNodeStats (
  ComponentSet *cs,
  Tree *tree,
  const Tree::Node *n,
  const DataComponentMap *dcm,
  const unsigned int w,
  const unsigned int h,
  float *ccimg,
  const float *refimg)
{
  ComponentSet::Component *cc, *cp;
  const unsigned int s = w*h;
  const Tree::Node *lf;
  unsigned int i, ii, cid, cnt;
  const float *por;
  float or[2], of[2], av[3];

  static unsigned int cccnt;
  // build the data for this component
  fprintf (stdout, "\rComponent %i/%i\0", cccnt++, cs->length ());

  // get the texture for this component
  if ((cc = cs->find (n->getCID ())) == NULL)
    return;

  memset (ccimg, 0, w*h*4*sizeof (float));

  // map the image
  for (i=0, ii=0; i<s; i++, ii+=4)
  {
    // copy the colour data
    ccimg[ii+0] = refimg[ii+0];
    ccimg[ii+1] = refimg[ii+1];
    ccimg[ii+2] = refimg[ii+2];
    ccimg[ii+3] = 0;

    // check if this pixel is a child of the current component
    cid = dcm->getDataComponent (i);

    if ((lf = tree->find (cid)) == NULL)
      continue;

    if (!lf->isAncestor (n))
      continue;

    // display this pixel in the image
    ccimg[ii+3] = 1.0f;
  }

  // upload the texture if this node is a leaf
#if 1
  if (n->isLeaf ())
    cc->setTextureData (w, h, ccimg, GL_RGBA, GL_FLOAT);
#else
//  char  buf[MAX_PATH] = {0};

//  sprintf_s (buf, MAX_PATH, "%i.png\0", n->getCID ());
//  PNGIO::writePixelBuffer (buf, ccimg, w, h, 8, 4);
#endif
  // otherwise, nodes don't have textures

  // build the outline (hull of this node)

  // build the origin of this node (center of gravity)
  // use the image data, not the stored texture so nodes that do not have
  // textures still have valid stats
  // init
  or[0]=0.0f;  or[1]=0.0f;
  av[0]=0.0f;  av[1]=0.0f;  av[2]=0.0f;
  for (i=0, ii=0, cnt=0; i<s; i++, ii+=4)
  {
    if (ccimg[ii+3] != 0.0f)
    {
      // component origin
      or[0] += (float)(i%w)/(float)w;
      or[1] += (float)(i/w)/(float)w;
      ++cnt;

      // component colour
      av[0] += ccimg[ii+0];
      av[1] += ccimg[ii+1];
      av[2] += ccimg[ii+2];
    }
  }

  or[0] /= (float)cnt;
  or[1] /= (float)cnt;
  av[0] /= (float)cnt;
  av[1] /= (float)cnt;
  av[2] /= (float)cnt;

  cc->setOrigin (or, 2);
  cc->setAverage (av, 3);

  // build the unit warp of this component
//  buildUnitWarp (w, h, ccimg, w, h, cc);

  // get the offset as the direction between this origin and the parent origin
  if (n->parent () == NULL)
  {
    of[0] = 0.0f;
    of[1] = 0.0f;
  }
  else
  {
    if ((cp = cs->find (n->parent ()->getCID ())) == NULL)
      return;

    por = cp->getOrigin ();

    of[0] = or[0]-por[0];
    of[1] = or[1]-por[1];
  }

  // normalise
//  l = sqrt (of[0]*of[0]+of[1]*of[1]);
//  of[0] /= l;
//  of[1] /= l;

  cc->setOffset (of, 2);

  // process the children
  if (n->left () != NULL)
  {
    getNodeStats (cs, tree, n->left (), dcm, w, h, ccimg, refimg);
  }
  if (n->right () != NULL)
  {
    getNodeStats (cs, tree, n->right (), dcm, w, h, ccimg, refimg);
  }
}

void buildComponentOffsets (
  ComponentSet *cs,
  Tree *tree,
  const DataComponentMap *dcm,
  const unsigned int refTexture)
{
  const Tree::Node *n;
  unsigned int t;
  float *ccimg, *rfimg;  // child buffer and reference image
  int w, h, s;

  // alloc
  glBindTexture (GL_TEXTURE_2D, refTexture);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);
  s = w*h;

  ccimg = new float[s*4];
  rfimg = new float[s*4];

  glGetTexImage (GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, rfimg);

  fprintf (stdout, "Building ccOffsets...\0");
  t = clock ();
#if 1
  n = tree->root ();
  getNodeStats (cs, tree, n, dcm, w, h, ccimg, rfimg);
#else
  // build the stuff
  for (i=0, cnt=0; i<tree->length (); i++)
  {
    if ((n = tree->get (i)) == NULL)
      continue;

    if ((p = n->parent ()) == NULL)
      continue;

    cc = cs->find (n->getCID ());
    cp = cs->find (p->getCID ());

    if ((cc == NULL) || (!cc->hasTexture ()) || (cp == NULL) || (!cp->hasTexture ()))
      continue;

    fprintf (stdout, "\rBuilding ccOffsets [%i/%i]...\0", cnt++, cs->length ());

    // get the image
    glBindTexture (GL_TEXTURE_2D, cc->getTexID ());
    glGetTexImage (GL_TEXTURE_2D, 0, GL_ALPHA, GL_FLOAT, ccimg);

    // get the parent image
    glBindTexture (GL_TEXTURE_2D, cp->getTexID ());
    glGetTexImage (GL_TEXTURE_2D, 0, GL_ALPHA, GL_FLOAT, cpimg);

    // count
    for (j=0, cccnt=0.0f, cpcnt=0.0f; j<s; j++)
    {
      if (ccimg[j] == 1.0f)
        cccnt+=1.0f;

      if (cpimg[j] == 1.0f)
        cpcnt+=1.0f;
    }

    // map the directions
    dm[0] = 0.0f;
    dm[1] = 0.0f;

    for (j=0; j<s; j++)
    {
      if (ccimg[j] != 1.0f)
      {
        continue;
      }

      p1[0] = (float)(j%w)/(float)s;
      p1[1] = (float)(j/w)/(float)s;
#if 1
      ds[0] = 0.0f;
      ds[1] = 0.0f;

      for (k=0; k<s; k++)
      {
        if (cpimg[j] != 1.0f)
        {
          continue;
        }

        p2[0] = (float)(k%w)/(float)s;
        p2[1] = (float)(k/w)/(float)s;

        ds[0] += (1.0f/cpcnt)*(p1[0]-p2[0]);
        ds[1] += (1.0f/cpcnt)*(p1[1]-p2[1]);
      }

      dm[0] += (1.0f/cccnt)*(ds[0]);
      dm[1] += (1.0f/cccnt)*(ds[1]);
#else
    //  dm[0] += (1.0f/cc)*(p1[0]);
    //  dm[1] += (1.0f/cc)*(p1[1]);
      dm[0] += p1[0];
      dm[1] += p1[1];
#endif
    }

    // update the offset
    l = sqrt (dm[0]*dm[0]+dm[1]*dm[1]);

    dm[0] /= l;
    dm[1] /= l;

    cc->setOffset (dm, 2);
  }
#endif

  fprintf (stdout, "\rBuilt ccOffsets in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
  // cleanup
  delete[] ccimg;
  delete[] rfimg;

  // build the tabs
  fprintf (stdout, "Building ccTabs...\0");
  t = clock ();

  getNodeTabs (cs, tree, tree->root (), 32, 32);

  fprintf (stdout, "\rBuild ccTab images in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
}

/*===============================================
 Images
===============================================*/
void buildDiscontinuityMap (
  DataComponentMap *dcm,
  Tree *tree,
  const unsigned int refTexture,
  const char *stub)
{
  int x, y, lbl[3];
  int w, h;
  unsigned char *buf;

  glBindTexture (GL_TEXTURE_2D, refTexture);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

  buf = new unsigned char[w*h*3];
  glGetTexImage (GL_TEXTURE_2D, 0, GL_RGB, GL_UNSIGNED_BYTE, buf);

#if 1
  // write the component edges to the image
  for (y=1; y<h; y++)
  {
    for (x=1; x<w; x++)
    {
      lbl[0] = dcm->getDataComponent ((y*w)+x);
      lbl[1] = dcm->getDataComponent (((y-1)*w)+x);
      lbl[2] = dcm->getDataComponent ((y*w)+(x-1));

      if ((lbl[0] != lbl[1]) || (lbl[0] != lbl[2]))
      {
        buf[(y*w*3)+(x*3)+0] = 255;
        buf[(y*w*3)+(x*3)+1] = 255;
        buf[(y*w*3)+(x*3)+2] = 255;
      }
    }
  }

#else
  unsigned int  cnt = dcm->uniques ();


#endif

  // write the image
  char  fname[MAX_PATH] = {0};

  sprintf_s (fname, MAX_PATH, "%s.png\0", stub);
  PNGIO::writePixelBuffer (fname, buf, w, h); 

  // cleanup
  delete[] buf;
}

/*===============================================
 Main functions
===============================================*/
#define FILE_OPEN              1
#define FILE_CLOSE             2
#define FILE_EXIT              3
#define FILE_SET_CACHE_DIR     4
#define FILE_LOAD_CACHE        5
#define FILE_GEN_CACHE         6
#define EDIT_SEL_ALL          10
#define EDIT_SEL_NONE         11
#define EDIT_COLLAPSE         12
#define BUILD_POINTSET        20
#define CLEAN_POINTSET        21
#define BUILD_COMPLETE_GRAPH  22
#define CLEAN_COMPLETE_GRAPH  23
#define BUILD_MST             24
#define CLEAN_MST             25
#define BALANCE_START         30
#define BALANCE_ITER          31
#define BALANCE_END           32
#define BALANCE_FAST          33
//#define PRUNE_NODES         24
#define LAYOUT_TRADITIONAL    40
#define LAYOUT_RADIAL         41
#define DRAWSTYLE_NODELINK    42
#define DRAWSTYLE_BLOCK       43
#define DRAWSTYLE_BLOCKTEX    44
#define PLOT_HIERARCHY        50
#define DRAW_BASIC            51
#define DRAW_INDENT           52
#define DRAW_REINGOLD         53
#define DRAW_GRUNDO           54
#define DRAW_AMDAM            55

#define APP_NAME "MST Clustering\0"

void createMainMenu (HWND hwnd)
{
  HMENU  root = CreateMenu ();
  HMENU  file = CreatePopupMenu ();
  HMENU  edit = CreatePopupMenu ();
  HMENU  opts = CreatePopupMenu ();

  AppendMenu (root, MF_POPUP|MF_STRING, (UINT_PTR)file, "&File\0");
  AppendMenu (file, MF_STRING, FILE_OPEN, "&Open...\0");
  AppendMenu (file, MF_STRING, FILE_CLOSE, "&Close\0");
  AppendMenu (file, MF_SEPARATOR, NULL, NULL);
  AppendMenu (file, MF_STRING, FILE_EXIT, "E&xit\0");

  AppendMenu (root, MF_POPUP|MF_STRING, (UINT_PTR)edit, "&Edit\0");
  AppendMenu (edit, MF_STRING, EDIT_SEL_ALL, "Select all\0");
  AppendMenu (edit, MF_STRING, EDIT_SEL_NONE, "Select none\0");
  AppendMenu (edit, MF_SEPARATOR, NULL, NULL);
  AppendMenu (edit, MF_STRING, EDIT_COLLAPSE, "Collapse node\0");

  AppendMenu (root, MF_POPUP|MF_STRING, (UINT_PTR)opts, "&Options\0");
  AppendMenu (opts, MF_STRING, BALANCE_START, "Balance tree\0");
  AppendMenu (opts, MF_SEPARATOR, 0, NULL);
  AppendMenu (opts, MF_STRING, LAYOUT_TRADITIONAL, "Traditional layout\0");
  AppendMenu (opts, MF_STRING, LAYOUT_RADIAL, "Radial Layout\0");
  AppendMenu (opts, MF_SEPARATOR, 0, NULL);
  AppendMenu (opts, MF_STRING, DRAWSTYLE_NODELINK, "Node-link style\0");
  AppendMenu (opts, MF_STRING, DRAWSTYLE_BLOCK, "Block style\0");
  AppendMenu (opts, MF_STRING, DRAWSTYLE_BLOCKTEX, "Tabsums style\0");
  AppendMenu (opts, MF_SEPARATOR, 0, NULL);
  AppendMenu (opts, MF_STRING, DRAW_BASIC, "Draw basic\0");
  AppendMenu (opts, MF_STRING, DRAW_INDENT, "Draw indented\0");
  AppendMenu (opts, MF_STRING, DRAW_REINGOLD, "Draw Reingold\0");
  AppendMenu (opts, MF_STRING, DRAW_GRUNDO, "Draw Grundo\0");
  AppendMenu (opts, MF_STRING, DRAW_AMDAM, "Draw Amsterdam\0");
  AppendMenu (opts, MF_SEPARATOR, 0, NULL);
  AppendMenu (opts, MF_STRING, PLOT_HIERARCHY, "Replot tree\0");

  SetMenu (hwnd, root);
}

LRESULT CALLBACK MainProc (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
  static unsigned int texid;
  static GLUtils::Camera camera;
  static ComponentSet *componentSet;
  static SegmentationTree *segTree;
  static AVLBalance *avb;
  static TreePlotter::PlotMethod plotMethod;
  static char cacheDir[MAX_PATH];

  switch (message)
  {
    case WM_CREATE:      // creation
    {
      createMainMenu (hwnd);

      if (!GLUtils::setupGLContext (hwnd))
        WinUtils::winErrorBox ("OpenGL context could not be created.\0");

      GLUtils::Text::createFonts (GetDC (hwnd), 10, "Verdana\0");

      // setup the camera
      float  defaultViewAngle[4] = {27.5, 27.5, 0.0, -3.0};
      float  defaultViewPos[3] = {0.0f, 0.0f, 0.0f};
      camera.setDefaultViewAngle (defaultViewAngle);
      camera.setDefaultViewPos (defaultViewPos);
      camera.reset ();

      plotMethod = TreePlotter::AMDAM;

      PostMessage (hwnd, WM_COMMAND, FILE_OPEN, 0);
      break;
    }
    case WM_MOVE:
    case WM_SIZE:      // Resizing
    {
      RECT  r;

      GetClientRect (hwnd, &r);

      glViewport (0, 0, r.right - r.left, r.bottom - r.top);
      break;
    }
    case WM_TIMER:      // intentional fall through to wm_paint
    {
      switch (wParam)
      {
      case 1:
        SendMessage (hwnd, WM_COMMAND, BALANCE_ITER, 0);
        break;
      case 2:
        segTree->buildHierarchyPlot (plotMethod);
        break;
      }
    }
    case WM_PAINT:      // Drawing (WM_TIMER and WM_PAINT are combined)
    {
      RECT  r;
      float  vp[4];

      GetClientRect (hwnd, &r);
      glGetFloatv (GL_VIEWPORT, vp);

      // openGL drawing
      glMatrixMode (GL_PROJECTION);  glLoadIdentity ();  glOrtho (0, r.right - r.left, r.bottom - r.top, 0, 0, 1);
      glMatrixMode (GL_MODELVIEW);  glLoadIdentity ();

      // clear
      glClearColor (1,1,0,0);
      glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

      // draw the texture
    //  if (segTree != NULL)
    //    segTree->drawTexture ();

      // draw the complete graph
    //  if (completeEdges != NULL)
    //    drawCompleteEdges (ps, psc, camera, completeEdges, completeEdgeCnt);

      // draw the 3d view of the points on the left hand side of the viewport
      if (segTree !=  NULL)
      {
        glPushAttrib (GL_VIEWPORT);

        glViewport (0, 0, vp[2]/2.0f, vp[3]);

        // draw the MST
        segTree->setupView (camera);
      //  segTree->drawMSTEdges ();
      //  segTree->drawPointSet (true);
      //  segTree->drawLeafComponents ();
        segTree->drawPoints ();
        segTree->cleanupView ();

        glPopAttrib ();
      }

      // draw the 2d tree on the right hand side of the viewport
      if (segTree != NULL)
      {
        glPushAttrib (GL_VIEWPORT);
        glViewport (vp[2]/2.0f, 0, vp[2]/2.0f, vp[3]);

        segTree->drawMSTTree ();

        glPopAttrib ();
      }

      // swap buffers
      SwapBuffers (GetDC (hwnd));

      // save a video
      if (avb != NULL)
      {
        const char  *dir  = "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\vids\\anim\\\0";
        const char  *fname = "avlbal\0";
        static unsigned int  frid;
        char  buf[MAX_PATH];

        sprintf_s (buf, MAX_PATH, "%s%s_%i.png\0", dir, fname, frid++);

        GLUtils::saveFrameBuffer (buf);
      }

      // end
      break;
    }
    case WM_COMMAND:    // Input commands
    {
      switch (LOWORD(wParam))
      {
        case FILE_OPEN:
        {
          char  fname[MAX_PATH] = {0};
          int    w, h;
#if 1
          char  fpath[MAX_PATH];    // buffer for file name

          if (!WinUtils::getOpenFileName (hwnd, fpath, MAX_PATH, "*.png\0"))
            break;
#else
          char  *fpath = TEST_IMAGE;
#endif
          glGenTextures (1, &texid);

          if (!GLUtils::readImageToTexture (fpath, texid, GL_RGB))
          {
            WinUtils::winErrorBox ("Unable to open file\0");
            glDeleteTextures (1, &texid);
            break;
          }

          // get the file name
          {
            char  *s, *e;

            s = strrchr (fpath, '\\');
            e = strrchr (fpath, '.');

            if (s == NULL)
              s = fpath;
            if (e == NULL)
              e = fpath + strlen (fpath);

            memcpy (fname, s+1, e-(s+1));
          }

          glBindTexture (GL_TEXTURE_2D, texid);
          glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
          glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

          fprintf (stdout, "SUC: Read file '%s' [%i x %i]\n\0", fname, w, h);

          // set the cache dir
          SendMessage (hwnd, WM_COMMAND, FILE_SET_CACHE_DIR, (LPARAM)fpath);

          // create the tree
          segTree = new SegmentationTree (texid);

          if (SendMessage (hwnd, WM_COMMAND, FILE_LOAD_CACHE, 0) == 0)
          {
            fprintf (stdout, "No cache found in '%s'\n\0", cacheDir);
            SendMessage (hwnd, WM_COMMAND, FILE_GEN_CACHE, 0);

            fprintf (stdout, "Loading newly cached data from '%s'\n\0", cacheDir);
            SendMessage (hwnd, WM_COMMAND, FILE_LOAD_CACHE, 0);
          }

          SendMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);

          // repaint
          PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case FILE_CLOSE:
        {
          if (glIsTexture (texid))
            glDeleteTextures (1, &texid);

          if (segTree != NULL)
          {
            delete segTree;
            segTree = NULL;
          }

          // redraw
          PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case FILE_EXIT:
        {
          SendMessage (hwnd, WM_CLOSE, 0, 0);
          break;
        }
        case FILE_SET_CACHE_DIR:
        {
          char  *path, *s, *e;
          char  fname[MAX_PATH] = {0};

          if (lParam == NULL)
            break;

          path = (char *)lParam;

          s = strrchr (path, '\\');
          e = strrchr (path, '.');

          if (s == NULL)
            s = path;
          if (e == NULL)
            e = path + strlen (path);

          memcpy (fname, s+1, e-(s+1));

          sprintf_s (cacheDir, MAX_PATH, "%s\\%s\\\0", CACHE_DIR, fname);

          // make the folder
          _mkdir (cacheDir);

          // make the ccs images folder
          char  buf[MAX_PATH];
          sprintf_s (buf, MAX_PATH, "%s%s\\\0", cacheDir, COMPONENTSET_IMG_SBDIR);

          _mkdir (buf);

          break;
        }
        case FILE_LOAD_CACHE:
        {
          DataComponentMap  *dcm;
          ComponentSet    *cs;
          Tree        *tree;

          // check if we can load the tree
          if ((tree = Tree::read (cacheDir, "t_prun\0")) == NULL)
          {
            return 0;
          }

          if ((dcm = DataComponentMap::read (cacheDir, "t_prun\0")) == NULL)
          {
            delete tree;
            return 0;
          }
#if 1
          // try to read the component set
          if ((cs = ComponentSet::read (cacheDir, "c_set\0")) == NULL)
          {
            // build components
            cs = buildComponents (cacheDir, texid, tree, dcm);
            buildComponentOffsets (cs, tree, dcm, texid);

            cs->write (cacheDir, "c_set\0");

            delete cs;

            cs = ComponentSet::read (cacheDir, "c_set\0");
          //  return 0;
          }
#else
          cs = NULL;
#endif
          segTree->setMSTHierarchy (tree);
          segTree->setDataComponentMap (dcm);
          segTree->setComponentSet (cs);

      //    buildDiscontinuityMap (dcm, tree, texid);
          return 1;
        }
        case FILE_GEN_CACHE:
        {
          DataComponentMap  *dcm;
          Tree        *tree;
          unsigned int    t, i, lc;
          TreePruner      *tp;

          // check if we only need to prune the tree (a full tree is cached)
          if ((tree = Tree::read (cacheDir, "t_full\0")) == NULL)
          {
            // construct the point set
            tree = (Tree*)SendMessage (hwnd, WM_COMMAND, BUILD_MST, 0);
          //  tree->buildWeights ();

            // write the complete tree to disk
            tree->write (cacheDir, "t_full\0");
          }

          if ((dcm = DataComponentMap::read (cacheDir, "t_full\0")) == NULL)
          {
            delete tree;
            break;
          }
#ifndef NO_PRUNE
          // prune the nodes
          fprintf (stdout, "Pruning tree...\0");
          t = clock ();
        //  segTree->getMSTHierarchy ()->pruneInternalNodes (256, segTree->getDataCIDMap (), segTree->getDataCIDMapLength (), true);
        //  tp = new TreePruner_Linear (segTree->getMSTHierarchy (), true);
          tp = new TreePruner_DFS (tree, false);

          tp->init ();
          tp->prune (PRUNE_DEPTH, dcm);

          delete tp;

          fprintf (stdout, "\rPruned tree in %0.2fs\n\0", (float)(clock()-t)/CLOCKS_PER_SEC);
#endif
          fprintf (stdout, "MST (CC=%i, RS=%i)\n\0", tree->length (), tree->root ()->size ());

          // write the tree
          tree->write (cacheDir, "t_prun\0");
          delete tree;

          // write the dcm
          dcm->write (cacheDir, "t_prun\0");
        //  buildDiscontinuityMap (dcm, tree, texid, "dcm_prun\0");
          delete dcm;

          // balance the tree
        //  fprintf (stdout, "Balancing tree...\0");
        //  t = clock ();

        //  tree = Tree::read (cacheDir, "t_prun\0");

        //  avb = new AVLBalance (tree);

        //  while (avb->next ())
        //    ;;

        //  delete avb;

        //  tree->write (cacheDir, "t_prun\0");
        //  delete tree;
        //  fprintf (stdout, "\rBalanced tree in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);

          break;
        }
        case EDIT_SEL_ALL:
        {
          if (segTree == NULL)
            break;

          segTree->setSelectedComponent (segTree->getMSTHierarchy ()->root ()->getCID ());
          break;
        }
        case EDIT_SEL_NONE:
        {
          if (segTree == NULL)
            break;

          segTree->setSelectedComponent (0xFFFFFFFF);
          break;
        }
        case EDIT_COLLAPSE:
        {
          if (segTree == NULL)
            break;

          segTree->collapseSelected ();
          break;
        }
        case BUILD_MST:
        {
          unsigned int  completeEdgeCount;
          unsigned int  *completeEdgeList;
          unsigned int  pointCount, pointDims;
          float      *pointSet;
          float      *weights;
          int        w, h;
          int        i, t;
          Tree      *tree = NULL;
          MSTree      *bruv = NULL;
          DataComponentMap  *dcm = NULL;

          if (segTree == NULL)
            break;

          glBindTexture (GL_TEXTURE_2D, texid);
          glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
          glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

          dcm = new DataComponentMap (w*h);
          constructPointSet (texid, &pointSet, pointCount, pointDims, false, dcm);
          dcm->write (cacheDir, "t_full\0");

        //  buildDiscontinuityMap (dcm, tree, texid, "dcm_full\0");

        //  segTree->setPointSet (pointSet, pointCount, pointDims);
#if 0
        //  completeEdgeCount = GraphUtils::completeGraphEdgeCount (pointCount);
          completeEdgeCount = GraphUtils::pixelNeighbourGraphEdgeCount (w, h);
          completeEdgeCount = GraphUtils::nearestNeighbourGraphEdgeCount (pointCount, 3);
          completeEdgeList = new unsigned int[completeEdgeCount*2];
          memset (completeEdgeList, 0, completeEdgeCount*2*sizeof (unsigned int));

          weights = new float[completeEdgeCount];
          memset (weights, 0, completeEdgeCount*sizeof (float));

        //  GraphUtils::buildCompleteGraph<unsigned int> (pointCount, completeEdgeList, completeEdgeCount);
        //  GraphUtils::buildPixelNeighbourGraph<unsigned int> (pointCount, w, h, completeEdgeList, completeEdgeCount);
          GraphUtils::buildNearestNeighbourGraph (pointSet, pointCount, completeEdgeList, completeEdgeCount, vertexWeight, 3);
          GraphUtils::fillEdgeWeights<float> (pointSet, pointCount, completeEdgeList, completeEdgeCount, weights, vertexWeight);

          fprintf (stdout, "Graph (V=%i, E=%i)\n\0", w*h, completeEdgeCount);

          // build a minimum spanning tree
          bruv = new MSTBoruvka (&tree, pointCount, completeEdgeList, weights, completeEdgeCount, true);

          fprintf (stdout, "MST[Bruv]...\0");
          t = clock ();

          for (bruv->init (), i=0; bruv->next (); i++)
            ;

          fprintf (stdout, "\rMST[Bruv] in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
          fprintf (stdout, "MST (CC=%i)\n\0", tree->length ());
        //  tree = MSTree::buildMST_Boruvka (segTree->getNumVertices (), completeEdgeList, weights, completeEdgeCount);

          delete[] weights;
#else
          bruv = new NNTree<float> (&tree, pointCount, true, vertexWeightd, pointSet);

          fprintf (stdout, "MST[NN]...\0");
          t = clock ();

          for (bruv->init (), i=0; bruv->next (); i++)
            ;

          fprintf (stdout, "\rMST[NN] in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
          fprintf (stdout, "MST (CC=%i)\n\0", tree->length ());
#endif

        //  delete tree;
        //  segTree->setMSTHierarchy (tree);

          // write the full tree
        //  SendMessage (hwnd, WM_COMMAND, FILE_WRITE_TREE, (LPARAM)"t_full\0");

          // FIXME: write the component data
          // cleanup
          delete[] pointSet;
          delete bruv;
          delete dcm;

        //  PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          return (LRESULT)tree;
        }
        case BALANCE_START:
        {
          if (segTree != NULL)
          {
            if (avb == NULL)
            {
              avb = new AVLBalance (segTree->getMSTHierarchy ());

              SetTimer (hwnd, 1, 10, NULL);
              SetTimer (hwnd, 2, 1000, NULL);
            }
          }
          break;
        }
        case BALANCE_ITER:
        {
          if ((segTree == NULL) || (avb == NULL))
            break;

          if (!avb->next ())
          {
            SendMessage (hwnd, WM_COMMAND, BALANCE_END, 0);
          }
          else
          {
            segTree->setSelectedComponent (avb->getCN ()->getCID ());
          }

        //  segTree->buildHierarchyPlot ();
        //  PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case BALANCE_END:
        {
          if (avb != NULL)
          {
            KillTimer (hwnd, 1);
            delete avb;
            avb = NULL;
          }
          if (segTree != NULL)
          {
            segTree->getMSTHierarchy ()->buildWeights ();
            segTree->buildHierarchyPlot (plotMethod);
          }
          break;
        }
        case BALANCE_FAST:
        {
          if (segTree != NULL)
          {
            if (avb == NULL)
            {
              avb = new AVLBalance (segTree->getMSTHierarchy ());

              while (avb->next ())
                ;;

              SendMessage (hwnd, WM_COMMAND, BALANCE_END, 0);
            }
          }
          break;
        }
        case PLOT_HIERARCHY:
        {
          unsigned int  t;

          if (segTree == NULL)
            break;

          fprintf (stdout, "Plotting tree...\0");
          t = clock ();
          segTree->buildHierarchyPlot (plotMethod);
          fprintf (stdout, "\rPlotted tree (%0.2fs)\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);

          PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case LAYOUT_TRADITIONAL:
        {
          if (segTree == NULL)
            return 0;

          segTree->setLayoutStyle (SegmentationTree::TRAD);
          PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          break;
        }
        case LAYOUT_RADIAL:
        {
          if (segTree == NULL)
            return 0;

          segTree->setLayoutStyle (SegmentationTree::RADIAL);
          PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          break;
        }
        case DRAWSTYLE_NODELINK:
        {
          if (segTree != NULL)
            segTree->setDrawStyle (SegmentationTree::NODELINK);

          break;
        }
        case DRAWSTYLE_BLOCK:
        {
          if (segTree != NULL)
            segTree->setDrawStyle (SegmentationTree::BLOCK);

          break;
        }
        case DRAWSTYLE_BLOCKTEX:
        {
          if (segTree != NULL)
            segTree->setDrawStyle (SegmentationTree::BLOCK_TEXTURES);

          break;
        }
        case DRAW_BASIC:
        {
          plotMethod = TreePlotter::DFS;
          PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          break;
        }
        case DRAW_INDENT:
        {
          plotMethod = TreePlotter::INDENT;
          PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          break;
        }
        case DRAW_REINGOLD:
        {
          plotMethod = TreePlotter::TR;
          PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          break;
        }
        case DRAW_GRUNDO:
        {
          plotMethod = TreePlotter::GRUNDO;
          PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          break;
        }
        case DRAW_AMDAM:
        {
          plotMethod = TreePlotter::AMDAM;
          PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          break;
        }
        default:
          break;
      }
      break;
    }
    case WM_KEYDOWN:
    {
      switch (wParam)
      {
        /*
        case VK_LEFT:
        case VK_RIGHT:
        {
          unsigned int  s;

          if (segTree == NULL)
            break;

          s = segTree->getSelectedComponent ();

          if ((wParam == VK_LEFT) && (s > 0))
            segTree->setSelectedComponent (s-1);
          else if ((wParam == VK_RIGHT) && (s<segTree->getMSTHierarchyCount ()))
            segTree->setSelectedComponent (s+1);

          break;
        }
        */
        case VK_UP:
        {
          unsigned int  s;

          if (segTree == NULL)
            break;

          s = segTree->getSelectedComponent ();

          segTree->selectParentComponent (s);
          break;
        }
        case VK_LEFT:
        {
          unsigned int  s;

          if (segTree == NULL)
            break;

          s = segTree->getSelectedComponent ();

          segTree->selectLeftComponent (s);
          break;
        }
        case VK_RIGHT:
        {
          unsigned int  s;

          if (segTree == NULL)
            break;

          s = segTree->getSelectedComponent ();

          segTree->selectRightComponent (s);
          break;
        }
        case VK_HOME:
        {
          PostMessage (hwnd, WM_COMMAND, EDIT_SEL_ALL, 0);
          break;
        }
        case VK_ESCAPE:
        {
          camera.reset ();
          PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case VK_SPACE:
        {
          PostMessage (hwnd, WM_COMMAND, EDIT_COLLAPSE, 0);
          break;
        }
        case 'B':
        case 'b':
        {
          PostMessage (hwnd, WM_COMMAND, BALANCE_START, 0);
          break;
        }
        case 'N':
        case 'n':
        {
          PostMessage (hwnd, WM_COMMAND, BALANCE_FAST, 0);
          break;
        }
        case VK_ADD:
        {
          float  s;

          if (segTree == NULL)
            break;

          s = segTree->getSeparation ();

          segTree->setSeparation (s*1.1f);
          break;
        }
        case VK_SUBTRACT:
        {
          float s;

          if (segTree == NULL)
            break;

          s = segTree->getSeparation ();

          segTree->setSeparation (s*0.9f);
          break;
        }
        default:
          break;
      }
      PostMessage (hwnd, WM_PAINT, 0, 0);
      return 0;
    }
    // close and exit
    case WM_DESTROY:
    {
      GLUtils::releaseGLContext (hwnd);
      return 0;
    }
    case WM_CLOSE:
    {
      SendMessage (hwnd, WM_COMMAND, FILE_CLOSE, 0);
      PostQuitMessage (0);
      return 0;
    }
    default:
    {
    //  if (camera.windowMessage (hwnd, message, wParam, lParam))
    //  {
    //    PostMessage (hwnd, WM_PAINT, 0, 0);
        break;
    //  }
    }
  }

  return DefWindowProc (hwnd, message, wParam, lParam);
}

int WINAPI WinMain (HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd)
{
  HWND    hwnd;
  WinUtils::Console    c;

  // setup the window
//  GLUtils::createGLWindowClass (hInstance, MainProc);
//  hwnd = GLUtils::createGLWindow (APP_NAME);
  WinUtils::createWindowClass (hInstance, MainProc, "CLSNME\0");
  hwnd = WinUtils::createWindow (APP_NAME, "CLSNME\0");

  WinUtils::showWindow (hwnd);

  WinUtils::resizeWindowClient (hwnd, 640, 480);

  WinUtils::messagePump (hwnd);
  WinUtils::destroyWindow (hwnd);

  return 1;
}
