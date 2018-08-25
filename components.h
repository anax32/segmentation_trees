#define  COMPONENTSET_FILE_TYPE  "ccs\0"
#define COMPONENTSET_IMG_SBDIR  "ccsimg\0"

class ComponentSet
{
public:
  class Component
  {
  public:
    static const unsigned int NO_TEX = 0xFFFFFFFF;    ///< indicates no texture for this component
    static const unsigned int MAX_DATA_DIMENSIONS = 3;  ///< maximum number of dimensions in the data-space
    static const unsigned int MAX_DATA_ATTRIBUTES = 4;  ///< maximum number of attributes for each sample

  protected:
    unsigned int cid;    ///< component id
    unsigned int texid;    ///< texture index
    unsigned int tabid;    ///< tab summary
    float offset[MAX_DATA_DIMENSIONS];  ///< direction to parent
    float origin[MAX_DATA_DIMENSIONS];  ///< center of mass for the component (in data space)
    float avg[MAX_DATA_ATTRIBUTES];    ///< average data value

  public:
    Component ()
      : cid (0), texid (NO_TEX), tabid (NO_TEX)
    {
      memset (offset, 0, MAX_DATA_DIMENSIONS*sizeof (float));
      memset (origin, 0, MAX_DATA_DIMENSIONS*sizeof (float));
      memset (avg, 0, MAX_DATA_DIMENSIONS*sizeof (float));
    }

    ~Component ()
    {}

    const bool hasTexture () const
    {
      return (texid != NO_TEX);
    }

    const bool hasTab () const
    {
      return (tabid != NO_TEX);
    }

    const unsigned int getCID () const
    {
      return cid;
    }

    const unsigned int getTexID () const
    {
      return texid;
    }

    const unsigned int getTabID () const
    {
      return tabid;
    }

    void setCID (const unsigned int newCID)
    {
      cid = newCID;
    }

    void setTextureData (
      const unsigned int w,
      const unsigned int h,
      const void *img,
      const unsigned int format = GL_RGBA,
      const unsigned int type = GL_FLOAT)
    {
      if (glIsTexture (texid) == GL_TRUE)
      {
        glDeleteTextures (1, &texid);
      }

      glGenTextures (1, &texid);
      glBindTexture (GL_TEXTURE_2D, texid);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, format, type, img);
    }

    void setTabData (
      const unsigned int w,
      const unsigned int h,
      const void *img,
      const unsigned int format = GL_RGBA,
      const unsigned int type = GL_FLOAT)
    {
      if (glIsTexture (tabid) == GL_TRUE)
      {
        glDeleteTextures (1, &tabid);
      }

      glGenTextures (1, &tabid);
      glBindTexture (GL_TEXTURE_2D, tabid);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      glTexImage2D (GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, format, type, img);
    }

    void setOffset (const float *offs, const unsigned int offsSize)
    {
      unsigned int  i;

      for (i=0; i<offsSize; i++)
      {
        offset[i] = offs[i];
      }
    }

    void setOrigin (const float *org, const unsigned int orgSize)
    {
      unsigned int  i;

      for (i=0; i<orgSize; i++)
      {
        origin[i] = org[i];
      }
    }

    void setAverage (const float *average, const unsigned int avSize)
    {
      unsigned int  i;

      for (i=0; i<avSize; i++)
      {
        avg[i] = average[i];
      }
    }

    const float *getOffset () const
    {
      return offset;
    }

    const float *getOrigin () const
    {
      return origin;
    }

    const float *getAverage () const
    {
      return avg;
    }

    // file io
    void read (FILE *f, const char *cacheDir)
    {
      unsigned int  w, h, ct, bd;
      char      buf[MAX_PATH];

      // read the component data
      fread (&cid, sizeof (unsigned int), 1, f);
      fread (offset, sizeof (float), MAX_DATA_DIMENSIONS, f);
      fread (origin, sizeof (float), MAX_DATA_DIMENSIONS, f);
      fread (avg, sizeof (float), MAX_DATA_ATTRIBUTES, f);

      // FIXME: this is a massive hack, for some mysterious reason y is inverted
      // between the read and write phase. Check it is not normalised somewhere.
      offset[1] = -offset[1];

      // try to read a texture from the cache
      sprintf (buf, "%s%s\\%i.png\0", cacheDir, COMPONENTSET_IMG_SBDIR, cid);

      if (PNGIO::readPixelBuffer (buf, NULL, w, h, bd, ct) != PNGIO::ERR_NOERR)
      {
        // no texture data
        texid = NO_TEX;
      }
      else if ((bd != 8) || (ct != PNGIO::RGBA))
      {
        fprintf (stdout, "ERR: Bad image format for component %i\n\0", cid);
        return;
      }
      else
      {
        unsigned char  *img;

        // alloc a buffer
        img = new unsigned char[w*h*4];
        memset (img, 0, w*h*4);

        PNGIO::readPixelBuffer (buf, img, w, h, bd, ct);

        setTextureData (w, h, img, GL_RGBA, GL_UNSIGNED_BYTE);

        delete[] img;
      }

      // try to read a tab from the cache
      sprintf (buf, "%s%s\\%i_tb.png\0", cacheDir, COMPONENTSET_IMG_SBDIR, cid);

      if (PNGIO::readPixelBuffer (buf, NULL, w, h, bd, ct) != PNGIO::ERR_NOERR)
      {
        tabid = NO_TEX;
      }
      else if ((bd != 8) || (ct != PNGIO::RGBA))
      {
        fprintf (stdout, "ERR: Bad image for tab %i\n\0", cid);
        return;
      }
      else
      {
        unsigned char  *img;

        // alloc a buffer
        img = new unsigned char[w*h*4];
        memset (img, 0, w*h*4);

        PNGIO::readPixelBuffer (buf, img, w, h, bd, ct);

        setTabData (w, h, img, GL_RGBA, GL_UNSIGNED_BYTE);

        delete[] img;
      }
    }

    void write (FILE *f, const char *cacheDir)
    {
      char      buf[MAX_PATH];

      // write the data
      fwrite (&cid, sizeof (unsigned int), 1, f);
      fwrite (offset, sizeof (float), MAX_DATA_DIMENSIONS, f);
      fwrite (origin, sizeof (float), MAX_DATA_DIMENSIONS, f);
      fwrite (avg, sizeof (float), MAX_DATA_ATTRIBUTES, f);

      // if we have a texture, write it to disk
      if (hasTexture ())
      {
        sprintf (buf, "%s%s\\%i.png\0", cacheDir, COMPONENTSET_IMG_SBDIR, cid);

        glBindTexture (GL_TEXTURE_2D, getTexID ());
        GLUtils::writeTextureToFile (buf, GL_RGBA);
      }

      // write the tab
      if (hasTab ())
      {
        sprintf (buf, "%s%s\\%i_tb.png\0", cacheDir, COMPONENTSET_IMG_SBDIR, cid);

        glBindTexture (GL_TEXTURE_2D, getTabID ());
        GLUtils::writeTextureToFile (buf, GL_RGBA);
      }
    }
  };

protected:
  const unsigned int  len;  ///< number of components
  const unsigned int  dsd;  ///< number of dimensions in the data-space (2 for images)

  Component      *cc;  ///< components

  void allocComponents ()
  {
    if (len > 0)
    {
      cc = new Component[len];
    }
  }
  void cleanComponents ()
  {
    if (cc != NULL)
    {
      delete[] cc;
      cc = NULL;
    }
  }
public:
  ComponentSet (const unsigned int numComponents, const unsigned int dataSpaceDimensions)
    : len (numComponents), dsd (dataSpaceDimensions), cc (NULL)
  {
    allocComponents ();
  }
  ~ComponentSet ()
  {
    cleanComponents ();
  }
  const unsigned int length () const
  {
    return len;
  }
  Component* get (const unsigned int index) const
  {
    return &cc[index];
  }
  Component* find (const unsigned int cid) const
  {
    unsigned int  i;

    for (i=0; i<length (); i++)
    {
      if (cc[i].getCID () == cid)
        return &cc[i];
    }

    return NULL;
  }
  // file io
  static ComponentSet* read (const char *cacheDir, const char *fname)
  {
    ComponentSet  *cs;
    unsigned int  i, _len, _dsd, t;
    FILE      *f;
    char      path[MAX_PATH];

    memset (path, 0, MAX_PATH);
    sprintf (path, "%s%s.ccs\0", cacheDir, fname);

    fprintf (stdout, "Reading components...\0");
    t = clock ();

    f = fopen (path, "rb\0");

    if (f == NULL)
    {
    //  fprintf (stdout, "ERR: File not found: '%s%s'\n\0", cacheDir, fname);
      return NULL;
    }

    fread (&_len, sizeof (unsigned int), 1, f);
    fread (&_dsd, sizeof (unsigned int), 1, f);

    cs = new ComponentSet (_len, _dsd);

    // read the components
    for (i=0; i<cs->length (); i++)
    {
      fprintf (stdout, "\rReading %i component\0", i);
      cs->get (i)->read (f, cacheDir);
    }

    fclose (f);

    fprintf (stdout, "\rRead %i components in %0.2fs\n\0", cs->length (), (float)(clock ()-t)/CLOCKS_PER_SEC);

    return cs;
  }
  void write (const char *cacheDir, const char *fname)
  {
    unsigned int  i, t;
    FILE      *f;
    char      path[MAX_PATH];

    memset (path, 0, MAX_PATH);
    sprintf (path, "%s%s.ccs\0", cacheDir, fname);

    fprintf (stdout, "Writing components...\0");
    t = clock ();

    f = fopen (path, "wb\0");

    if (f == NULL)
    {
      fprintf (stdout, "ERR: Unable to write to '%s'\n\0", fname);
      return;
    }

    // write the header
    fwrite (&len, sizeof (unsigned int), 1, f);
    fwrite (&dsd, sizeof (unsigned int), 1, f);

    // write the components
    for (i=0; i<length (); i++)
    {
      fprintf (stdout, "\rWriting component %i...\0", i);
      cc[i].write (f, cacheDir);
    }

    fclose (f);

    fprintf (stdout, "\rWrote %i components in %0.2fs\n\0", length (), (float)(clock ()-t)/CLOCKS_PER_SEC);
  }
};
