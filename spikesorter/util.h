
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#define DIE (fprintf (stderr, "fatal error in %s at line %d \n", __FILE__, __LINE__), die(), 0)

#define TMALLOC(buf, n) ((buf = malloc ((n) * sizeof *(buf)))								\
			 && (malloc_debug == 0 || fprintf (stderr, "%d %d %ld\n", __LINE__, 0, (long)buf))) || DIE

#define TCALLOC(buf, n) ((buf = calloc ((n), sizeof *(buf)))								\
			 && (malloc_debug == 0 || fprintf (stderr, "%d %d %ld\n", __LINE__, 0, (long)buf))) || DIE

#define TREALLOC(buf, n) ((malloc_debug == 0 || fprintf (stderr, "%d %ld", __LINE__, (long)buf))	\
			  && (buf = realloc (buf, (n) * sizeof *(buf)))			\
			  && (malloc_debug == 0 || fprintf (stderr, " %ld\n", (long)buf))) || DIE


extern int malloc_debug;
void die (void);
