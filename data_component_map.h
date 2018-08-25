class DataComponentMap
{
protected:
  unsigned int  len;
  unsigned int  *cid;

  void allocCIDMap ()
  {
    cid = new unsigned int[len];
    memset (cid, 0, len*sizeof (unsigned int));
  }

  void cleanCIDMap ()
  {
    if (cid != NULL)
    {
      delete[] cid;
      cid = NULL;
      len = 0;
    }
  }

public:
  DataComponentMap (const unsigned int numPoints)
    : len (numPoints), cid (NULL)
  {
    allocCIDMap ();
  }

  DataComponentMap ()
    : len (0), cid (NULL)
  {}

  virtual ~DataComponentMap ()
  {
    cleanCIDMap ();
  }

  void setDataComponent (const unsigned int did, const unsigned int ncid)
  {
    cid[did] = ncid;
  }

  const unsigned int getDataComponent (const unsigned int did) const
  {
    return cid[did];
  }

  const unsigned int length () const
  {
    return len;
  }

  const unsigned int uniques () const
  {
    unsigned int  i, j, cnt;

    cnt = 0;

    for (i=0; i<len; i++)
    {
      for (j=0; j<len; j++)
      {
        if (cid[i]==cid[j])
          break;
      }

      if (j==len)
        cnt++;
    }

    return cnt;
  }

  // io
  static DataComponentMap *read (const char *cacheDir, const char *fname)
  {
    DataComponentMap *dcm;
    FILE *f;
    char path[MAX_PATH];
    unsigned int l, t;

    memset (path, 0, MAX_PATH);
    sprintf (path, "%s%s.dcm\0", cacheDir, fname);

    fprintf (stdout, "Reading DataComponentMap...\0");
    t = clock ();

    f = fopen (path, "rb\0");

    if (f == NULL)
      return NULL;

    fread (&l, sizeof (unsigned int), 1, f);

    if (l == 0)
      return NULL;

    dcm = new DataComponentMap (l);

    fread (dcm->cid, sizeof (unsigned int), dcm->length (), f);

    fclose (f);

    fprintf (stdout, "\rRead DataComponentMap in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);

    return dcm;
  }

  void write (const char *cacheDir, const char *fname)
  {
    FILE *f;
    char path[MAX_PATH];
    unsigned int t;

    memset (path, 0, MAX_PATH);
    sprintf (path, "%s%s.dcm\0", cacheDir, fname);

    fprintf (stdout, "Writing DataComponentMap...\0");
    t = clock ();

    f = fopen (path, "wb\0");

    if (f == NULL)
      return;

    fwrite (&len, sizeof (unsigned int), 1, f);
    fwrite (cid, sizeof (unsigned int), len, f);
    fclose (f);

    fprintf (stdout, "\rWrote DataComponentMap in %0.2fs\n\0", (float)(clock ()-t)/CLOCKS_PER_SEC);
  }
};
