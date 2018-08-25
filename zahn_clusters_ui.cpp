#pragma comment(lib, "../../Image/PNG/libpng13.lib")  // libpng
//#pragma comment(lib, "../../glew/lib/glew32.lib")

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <assert.h>
#include <iostream>

#define OEMRESOURCE
#define NOMINMAX
#include <windows.h>
#include <windowsx.h>

#include "../../Utils/Patterns/Persistable.h"
#include "../../Utils/Structures/ParameterSet.h"
#include "../../Utils/Callbacks/Callback.h"
#include "../../Utils/Win/Utils.h"  // console
#include "../../Utils/Win/Messages.h"
#include "../../Utils/Win/Console.h"
#include "../../Image/pngreader.h"  // libpng
#include "../../Utils/ImageUtils.h"
#include "../../Utils/GraphUtils.h"

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

/*===============================================
 Construct a graph from a gl texture
===============================================*/
float vertexWeight (const float *vertices, const unsigned int v1, const unsigned int v2)
{
  float  x, y, d;

  x = abs (vertices[(v1*2)+0]-vertices[(v2*2)+0]);
  y = abs (vertices[(v1*2)+1]-vertices[(v2*2)+1]);

//  d = sqrt (x*x+y*y);
  d = pow (x, 5.0f) + pow (y, 5.0f);

  d = exp (x) + exp (y);

  return d;
}

void thresholdImage (unsigned int texture)
{
  unsigned int  i, j;
  unsigned int  w, h;
  unsigned char  *img, mn, mx, threshold;

  // first construct a point set of the texture
  glBindTexture (GL_TEXTURE_2D, texture);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, (GLint*)&w);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, (GLint*)&h);

  // allocate space
  img = new unsigned char[w*h*1];
  glGetTexImage (GL_TEXTURE_2D, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, img);

  // threshold the image
  mn = 0;
  mx = 255;
  threshold = 128;
  ImageUtils::threshold<unsigned char> (img, w, h, 1, &mn, &mx, &threshold);

  // remove any black pixels which have more than one black neighbour
  for (i=1; i<h-1; i++)
  {
    for (j=1; j<w-1; j++)
    {
      if (img[(i*w)+j] == 255)
        continue;

      // if all are black, leave it
      if ((img[(i*w)+(j-1)] == 0) && (img[(i*w)+(j+1)] == 0) &&
          (img[((i-1)*w)+j] == 0) && (img[((i+1)*w)+j] == 0) &&
          (img[((i-1)*w)+(j-1)] == 0) && (img[((i-1)*w)+(j+1)] == 0) &&
          (img[((i+1)*w)+(j-1)] == 0) && (img[((i+1)*w)+(j+1)] == 0))
      {
        continue;
      }

      // if some are not black, set this to white
      if ((img[(i*w)+(j-1)] == 0) || (img[(i*w)+(j+1)] == 0) ||
          (img[((i-1)*w)+j] == 0) || (img[((i+1)*w)+j] == 0) ||
          (img[((i-1)*w)+(j-1)] == 0) || (img[((i-1)*w)+(j+1)] == 0) ||
          (img[((i+1)*w)+(j-1)] == 0) || (img[((i+1)*w)+(j+1)] == 0))
      {
        img[(i*w)+j] = 255;
      }
    }
  }

  // reload the texture
  glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, img);

  delete[] img;
  img = NULL;
}

void constructPointSet (unsigned int texture, float **pointSet, unsigned int& pointCount)
{
  // build a point set from a thresholded texture
  unsigned int  i, j, c, k;
  unsigned int  w, h;
  unsigned char  *img;
  float      *p;

  // first construct a point set of the texture
  glBindTexture (GL_TEXTURE_2D, texture);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, (GLint*)&w);
  glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, (GLint*)&h);

  // allocate space
  img = new unsigned char[w*h];
  glGetTexImage (GL_TEXTURE_2D, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, img);

  // first count the number of points
  for (i=0, c=0; i<w*h; i++)
  {
    if (img[i] == 0)
      ++c;
  }

  // allocate space
  p = new float[c*2];
  memset (p, 0, c*2*sizeof (float));

  for (i=0, k=0; i<h && k<c; i++)
  {
    for (j=0; j<w && k<c; j++)
    {
      if (img[(i*w)+j] == 0)
      {
        p[(k<<1)+0] = ((float)j)-((float)(w>>1));
        p[(k<<1)+1] = ((float)h-i)-((float)(h>>1));
        ++k;
      }
    }
  }

  delete[] img;
  img = NULL;

  *pointSet = p;
  pointCount = k;
}

void buildAvgVertWeights (
  float *avgw,
  const unsigned int V,
  const unsigned int *edges,
  const float *w,
  const unsigned int E)
{
  unsigned int  i, s, d;
  unsigned int  *cnt;

  cnt = new unsigned int[V];
  memset (cnt, 0, V*sizeof (unsigned int));

  // sum
  for (i=0; i<E; i++)
  {
    s = edges[(i*2)+0];
    d = edges[(i*2)+1];

    avgw[s] += w[i];
    avgw[d] += w[i];
    cnt[s]++;
    cnt[d]++;
  }

  // avg
  for (i=0; i<V; i++)
  {
    avgw[i] /= (float)cnt[i];
  }

  // cleanup
  delete[] cnt;
}

/*===============================================
 Main functions
===============================================*/
#define FILE_OPEN            1
#define FILE_CLOSE           2
#define FILE_EXIT            3
#define IMG_THRESH           4
#define IMG_POINT_SET        5
#define IMG_CMPLT_GRAPH      6
#define IMG_MST_GRAPH        7
#define BALANCE_START       10
#define BALANCE_ITER        11
#define BALANCE_END         12
#define BALANCE_FAST        23
#define LAYOUT_TRADITIONAL  30
#define LAYOUT_RADIAL       31
#define PLOT_HIERARCHY      40
#define DRAW_BASIC          41
#define DRAW_INDENT         42
#define DRAW_REINGOLD       43

#define APP_CLASS "ZHNCLS\0"
#define APP_NAME "Gestalt Clusters (Zahn1971)\0"
#define TEST_IMAGE  "C:\\Users\\csed\\Documents\\Projects\\Gestalt\\data\\Zahn-two-sets-2.png\0"

void createMainMenu (HWND hwnd)
{
  HMENU  root = CreateMenu ();
  HMENU  file = CreatePopupMenu ();
  HMENU  image = CreatePopupMenu ();
  HMENU  opts = CreatePopupMenu ();

  AppendMenu (root, MF_POPUP|MF_STRING, (UINT_PTR)file, "&File\0");
  AppendMenu (file, MF_STRING, FILE_OPEN, "&Open...\0");
  AppendMenu (file, MF_STRING, FILE_CLOSE, "&Close\0");
  AppendMenu (file, MF_SEPARATOR, NULL, NULL);
  AppendMenu (file, MF_STRING, FILE_EXIT, "E&xit\0");

  AppendMenu (root, MF_POPUP|MF_STRING, (UINT_PTR)image, "&Image\0");
  AppendMenu (image, MF_STRING, IMG_THRESH, "&Threshold\0");
  AppendMenu (image, MF_STRING, IMG_POINT_SET, "&Point set\0");
  AppendMenu (image, MF_STRING, IMG_CMPLT_GRAPH, "&Complete graph\0");
  AppendMenu (image, MF_STRING, IMG_MST_GRAPH, "&MST\0");

  AppendMenu (root, MF_POPUP|MF_STRING, (UINT_PTR)opts, "Options\0");
  AppendMenu (opts, MF_STRING, LAYOUT_TRADITIONAL, "Traditional layout\0");
  AppendMenu (opts, MF_STRING, LAYOUT_RADIAL, "Radial Layout\0");
  AppendMenu (opts, MF_SEPARATOR, 0, NULL);
  AppendMenu (opts, MF_STRING, DRAW_BASIC, "Draw basic\0");
  AppendMenu (opts, MF_STRING, DRAW_INDENT, "Draw indented\0");
  AppendMenu (opts, MF_STRING, DRAW_REINGOLD, "Draw Reingold\0");

  SetMenu (hwnd, root);
}

LRESULT CALLBACK MainProc (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
  static unsigned int    texture, refTexture;
  static GLUtils::Camera  camera;
  static SegmentationTree  *segTree;
  static AVLBalance    *avb;
  static TreePlotter::PlotMethod  plotMethod;

  switch (message)
  {
    case WM_CREATE:      // creation
    {
      createMainMenu (hwnd);

      if (!GLUtils::setupGLContext (hwnd))
        MessageBox (hwnd, "OpenGL context could not be created.\0", "Error\0", MB_OK|MB_ICONERROR);

      GLUtils::Text::createFonts (GetDC (hwnd));

      // setup the camera
      float  defaultViewAngle[4] = {27.5, 27.5, 0.0, -3.0};
      float  defaultViewPos[3] = {0.0f, 0.0f, 0.0f};
      camera.setDefaultViewAngle (defaultViewAngle);
      camera.setDefaultViewPos (defaultViewPos);
      camera.reset ();

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
      glClearColor (0,0,0,0);
      glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

      // draw the texture
      if (segTree != NULL)
        segTree->drawTexture ();

      // draw the complete graph
    //  if (completeEdges != NULL)
    //    drawCompleteEdges (ps, psc, camera, completeEdges, completeEdgeCnt);

      // draw the 3d view of the points on the left hand side of the viewport
      if (segTree != NULL)
      {
        glPushAttrib (GL_VIEWPORT);

        glViewport (0, 0, vp[2]/2.0f, vp[3]);

        // draw the MST
        segTree->setupView (camera);
      //  segTree->drawMSTEdges ();
      //  segTree->drawPointSet (false);
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

      // end
      break;
    }
    case WM_COMMAND:    // Input commands
    {
      switch (LOWORD(wParam))
      {
        case FILE_OPEN:
        {
#if 1
          char fname[MAX_PATH];    // buffer for file name

          if (!WinUtils::getOpenFileName (hwnd, fname, MAX_PATH, "*.png\0"))
            break;
#else
          char *fname = TEST_IMAGE;
#endif
          // close the current file
          SendMessage (hwnd, WM_COMMAND, FILE_CLOSE, 0);

          // make an opengl texture
          glGenTextures (1, &texture);

          if (!GLUtils::readImageToTexture (fname, texture, GL_LUMINANCE))
          {
            glDeleteTextures (1, &texture);
            MessageBox (hwnd, "Could not get image details.\0", "Error\0", MB_OK|MB_ICONERROR);
            break;
          }
          else
          {
            int  w, h;
            glBindTexture (GL_TEXTURE_2D, texture);
            glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
            glGetTexLevelParameteriv (GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

            fprintf (stdout, "SUC: Opened '%s' [%i x %i x 1]\n\0", fname, w, h);
          }

          glGenTextures (1, &refTexture);

          if (!GLUtils::readImageToTexture (fname, refTexture, GL_LUMINANCE))
          {
            glDeleteTextures (1, &refTexture);
            MessageBox (hwnd, "Could not get image details.\0", "Error\0", MB_OK|MB_ICONERROR);
          }

          segTree = new SegmentationTree (refTexture);

          // process the image
          SendMessage (hwnd, WM_COMMAND, IMG_THRESH, 0);
          SendMessage (hwnd, WM_COMMAND, IMG_POINT_SET, 0);
          SendMessage (hwnd, WM_COMMAND, IMG_CMPLT_GRAPH, 0);
          SendMessage (hwnd, WM_COMMAND, IMG_MST_GRAPH, 0);

          // repaint
          PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case FILE_CLOSE:
        {
          if (glIsTexture (texture))  {glDeleteTextures (1, &texture);}
          if (glIsTexture (refTexture))  {glDeleteTextures (1, &refTexture);}

          if (segTree != NULL)  {delete segTree; segTree = NULL;}
          // redraw
          PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case FILE_EXIT:
        {
          SendMessage (hwnd, WM_CLOSE, 0, 0);
          break;
        }
        case IMG_THRESH:
        {
          if (!glIsTexture (texture))
          {
            MessageBox (hwnd, "A texture must be loaded.\0", "Error\0", MB_OK|MB_ICONERROR);
            break;
          }

          thresholdImage (texture);

          PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case IMG_POINT_SET:
        {
          if (!glIsTexture (texture))
          {
            MessageBox (hwnd, "A thresholded texture must be loaded.\0", "Error\0", MB_OK|MB_ICONERROR);
            break;
          }

          if (segTree == NULL)
            break;

          unsigned int  psc;
          float      *ps;

          constructPointSet (texture, &ps, psc);

          segTree->setPointSet (ps, psc, 2);

          delete[] ps;

          PostMessage (hwnd, WM_PAINT, 0, 0);
          break;
        }
        case IMG_CMPLT_GRAPH:
        {
          if (segTree == NULL)
            break;

          unsigned int  *e, ecnt;
          float      *w;

          // get the number of edges
        //  ecnt = GraphUtils::completeGraphEdgeCount (segTree->getNumVertices ());
          ecnt = GraphUtils::nearestNeighbourGraphEdgeCount (segTree->getNumVertices (), 5);

          // allocate space
          e = new unsigned int[ecnt*2];
          w = new float[ecnt];

          memset (e, 0, ecnt*2*sizeof (unsigned int));
          memset (w, 0, ecnt*sizeof (float));

          // construct
#if 1
        //  GraphUtils::buildCompleteGraph<unsigned int> (segTree->getNumVertices (), e, ecnt);
          GraphUtils::buildNearestNeighbourGraph<float> (
            segTree->getVertices (),
            segTree->getNumVertices (),
            e,
            ecnt,
            vertexWeight,
            5);
#else
          GraphUtils::buildPixelNeighbourGraph<unsigned int> (
            segTree->getNumVertices (),
            64,
            64,
            e,
            ecnt);
#endif
          // get the weights
          GraphUtils::fillEdgeWeights<float> (
            segTree->getVertices (),
            segTree->getNumVertices (),
            e,
            ecnt,
            w,
            vertexWeight);
#if 0
          // build the mst
          unsigned int *mstEdges, mstEdgeCnt;
          float *mstWeights;
          unsigned int *mstHierarchy, mstHierarchyCnt;

          mstEdgeCnt = 0;
          mstHierarchyCnt = 0;

          // construct
        //  constructMSTGraph (psc, completeEdges, completeWeights, completeEdgeCnt, &mstEdges, &mstWeights, mstEdgeCnt);
        //  constructMSTGraph_Kruskal (psc, completeEdges, completeWeights, completeEdgeCnt, &mstEdges, &mstWeights, mstEdgeCnt);
        //  GraphUtils::buildMST_Kruskals (segTree->getNumVertices (), e, w, ecnt, &mstEdges, &mstWeights, mstEdgeCnt, &mstHierarchy, mstHierarchyCnt);
          GraphUtils::buildMST_KruskalsAgg (
            segTree->getNumVertices (),
            e,
            w,
            ecnt,
            &mstEdges,
            &mstWeights,
            mstEdgeCnt,
            &mstHierarchy,
            mstHierarchyCnt);

          segTree->setMSTEdges (mstEdges, mstEdgeCnt, mstWeights);
          segTree->setMSTHierarchy (mstHierarchy, mstHierarchyCnt);

          // cleanup
          delete[] mstEdges;
          delete[] mstWeights;
          delete[] mstHierarchy;
#else
          int i, t;
          Tree *tree = NULL;
          MSTree *bruv = NULL;

          bruv = new MSTBoruvka (&tree, segTree->getNumVertices (), e, w, ecnt);

          fprintf (stdout, "MST[Bruv]...\n\0");
          t = clock ();

          bruv->init ();

          for (i=0; bruv->next (); i++)
          {
            fprintf (stdout, "MST[Bruv] %i [%i]\n\0", i, tree->countRoots ());
          }

          delete bruv;

          fprintf (stdout, "MST[Bruv] in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);

        //  tree = MSTree::buildMST_Boruvka (segTree->getNumVertices (), completeEdgeList, weights, completeEdgeCount);

        //  delete tree;
          segTree->setMSTHierarchy (tree);
#endif
          delete[] e;
          delete[] w;

          // redraw
          PostMessage (hwnd, WM_COMMAND, PLOT_HIERARCHY, 0);
          break;
        }
        case BALANCE_START:
        {
          if (segTree != NULL)
          {
            if (avb == NULL)
            {
              avb = new AVLBalance (segTree->getMSTHierarchy ());
              SetTimer (hwnd, 1, 10, NULL);
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

          segTree->buildHierarchyPlot (plotMethod);
          PostMessage (hwnd, WM_PAINT, 0, 0);
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
            segTree->buildHierarchyPlot (plotMethod);

          PostMessage (hwnd, WM_PAINT, 0, 0);
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
          if (segTree == NULL)
            break;

          segTree->setSelectedComponent (segTree->getMSTHierarchy ()->root ()->getCID ());
          break;
        }
        case VK_ESCAPE:
        {
          camera.reset ();
          PostMessage (hwnd, WM_PAINT, 0, 0);
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
        default:
          break;
      }
      PostMessage (hwnd, WM_PAINT,0,0);
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
    //    PostMessage (hwnd, WM_PAINT, 0, 0);
      break;
    }
  }

  return DefWindowProc (hwnd, message, wParam, lParam);
}

int WINAPI WinMain (HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd)
{
  HWND hwnd;
  WinUtils::Console  c;

  // setup the window
  WinUtils::createWindowClass (hInstance, MainProc, APP_CLASS);
  hwnd = WinUtils::createWindow (APP_NAME, APP_CLASS);

  WinUtils::showWindow (hwnd);
  WinUtils::messagePump (hwnd);
  WinUtils::destroyWindow (hwnd);

  return 1;
}
