/*
 * stringutils.c
 * functions to manipulate strings
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: stringutils.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/stringutils.c $
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <assert.h>
 #include "stringutils.h"
 #include "basic-types.h"


void
printSubseq(char* seq, Uint start, Uint end) {
  int i;
  assert(end <= strlen(seq));
  for(i=start; i <= end; i++) {
    printf("%c", seq[i]);
  }
}

char* strrev(char *str, Uint len){
  int end = len-1;
  int start = 0;

  while (start<end) {
    str[start] ^= str[end];
    str[end] ^= str[start];
    str[start] ^= str[end];
    start++;
    end--;
  }
  return str;
}

char* strtok_bl(char *s, char *delim, char **saveptr) {

  char *ret;

  /* init */
  if (s == NULL){
    s = *saveptr;
    if (s == NULL){
      *saveptr = NULL;
      return NULL;
    }
  }
  /* skip delims at begin */
  while(*s && strchr(delim, *s)){
    s++;
  }

  if (*s == 0){
    *saveptr = NULL;
    return NULL;
  }

  /* locate next delim or end */
  ret = s;
  while(*s && !strchr(delim, *s)){
    s++;
  }

  if (*s == 0){
    *saveptr = NULL;
  }
  else {
    *s = 0;
    s++;
    *saveptr = s;
  }  

  return ret;
}


stringset_t* tokensToStringset(void *space, char* delim, char* toTokens, 
                                 Uint len){
	Uint toklen;
	char* token;
        char* saveptr;
	char* buffer;
	stringset_t *set;
	
	set = ALLOCMEMORY(space, NULL, stringset_t, 1);
	set->noofstrings = 0;
	set->strings = NULL;
									
	if (toTokens == NULL || len == 0)
	return set;
													
	buffer = ALLOCMEMORY(space, NULL, char, len+1);
	buffer = memcpy(buffer, toTokens, len+1);
        buffer[len] = 0;
																
	if (buffer == NULL) {
		fprintf(stderr, "copy tokenstring %s to buffer failed.\n", toTokens);
		exit(-1);	
	}
	
	token = strtok_bl(buffer, delim, &saveptr);
	
	while(token != NULL) {
		
			toklen = strlen(token);
			set->noofstrings++;
			set->strings = ALLOCMEMORY(space, set->strings, string_t, set->noofstrings);
			set->strings[set->noofstrings-1].str = ALLOCMEMORY(space, NULL, char, toklen+1);
			set->strings[set->noofstrings-1].str = memcpy(set->strings[set->noofstrings-1].str, token, toklen);
			set->strings[set->noofstrings-1].str[toklen]='\0'; 
			set->strings[set->noofstrings-1].len = toklen;
			token = strtok_bl(NULL, delim, &saveptr);
	}

	FREEMEMORY(space, buffer);
	return set;
}


char* strtrimquote (void *spacetab, char *toTrim, Uint *len) {
	Uint i=0;
	int start=-1;
	int end =-2;	
	char* trimmed = NULL;
	
	for(i=0; i < *len; i++) {
	  
		if(ISQUOTE(toTrim[i])) {	
			continue;
		}
		else if(start == -1) {
			start=i;
			end=i;
		}
		else 
			end=i;
	}

	if(start >= 0) {
		trimmed = ALLOCMEMORY(spacetab, NULL, char, (end-start)+2);
   		memmove(trimmed, &toTrim[start], (end-start)+1);
		trimmed[(end-start)+1]='\0';
	}
	
	*len = (end-start)+1;
	return trimmed;

}


char* strtrim(void *spacetab, char* toTrim, Uint *len) {	
	Uint i=0;
	int start=-1;
	int end =-2;	
	char* trimmed = NULL;
	
	for(i=0; i < *len; i++) {
	  
		if(ISWHITESPACE(toTrim[i])) {	
			continue;
		}
		else if(start == -1) {
			start=i;
			end=i;
		}
		else 
			end=i;
	}

	if(start >= 0) {
		trimmed = ALLOCMEMORY(spacetab, NULL, char, (end-start)+2);
   		memmove(trimmed, &toTrim[start], (end-start)+1);
		trimmed[(end-start)+1]='\0';
	}
	
	*len = (end-start)+1;
	return trimmed;
}

char* strclip(void *spacetab, char* str, Uint *len){
  Uint i;
  char *tmp;
  // clip string after first whitespace
  for (i = 0; i < *len; i++){
    if (ISWHITESPACE(str[i])){
      break;
    }
  }
  
  tmp = ALLOCMEMORY(spacetab, NULL, char, i + 1);
  memmove(tmp, str, i);
  tmp[i] = '\0';
  *len = i;
  return tmp;
}

void strconvert(char *seq, Uint len, char orig, char replace){
  Uint i;
  for (i = 0; i < len; i++){
    if (seq[i] == orig){
      seq[i] = replace;
    }
  }
}

char* concat(void *spacetab, char* strA, char* strB, int lenA, int lenB) {

	if(strB == NULL || lenB == 0) 
	  		return strA;
	if(strA == NULL || lenA == 0)
		  		return strB;
			
	strA=ALLOCMEMORY(spacetab, strA, char, lenA+lenB+1);
	memmove(&strA[lenA], strB, lenB);
					
	strA[lenA+lenB]='\0';
					
	return strA;
}

char* concatdelim(void *spacetab, char* strA, char* strB, int lenA, int lenB, char delim) {

	if(strB == NULL || lenB == 0) 
	  		return strA;
	if(strA == NULL || lenA == 0)
		  	return strB;
			
	strA=ALLOCMEMORY(spacetab, strA, char, lenA+lenB+2);
	strA[lenA]=delim;
	memmove(&strA[lenA+1], strB, lenB);
						
	strA[lenA+lenB+1]='\0';
						
	return strA;						
}


void destructStringset(void *space, stringset_t *s) {
  Uint i;

  if (s->strings) {
    for(i=0; i < s->noofstrings; i++) {
      if(s->strings[i].str != NULL)	
        FREEMEMORY(space, s->strings[i].str);
    }

    FREEMEMORY(space, s->strings);
  }

  FREEMEMORY(space, s);
}


stringset_t *initStringset(void *space) {
	stringset_t *set;

	set = ALLOCMEMORY(space, NULL, stringset_t, 1);
	set->noofstrings=0;
	set->strings = NULL;

	return set;	
}

void addString(void *space, stringset_t *set, char *string, Uint len) {

	set->noofstrings++;
	set->strings=ALLOCMEMORY(space, set->strings, string_t, set->noofstrings);
	set->strings[set->noofstrings-1].str = string;
	set->strings[set->noofstrings-1].len = len;
}



/*---------------------------------- strrev ----------------------------------
 *    
 * reversing a string
 * 
 */
 
char*
strreverse(char *s, Uint len)
{	
    Uint i;
	char resc;
	
  	for(i=0; i < (len/2); i++) {
		resc = s[i];
		s[i] = s[len-1-i];
		s[len-1-i] = resc;
	}
	
	return s;
}

/* -------------------------------- my_itoa  ---------------------------------
 *    
 * just in case that there's no itoa
 * 
 */
 
char*
my_itoa (int value, char *buffer, Uint radix)
{
  	const char* base ="0123456789abcdef";
    int i=0;

	if (value == 0) {
	  buffer[0]=base[0];
	  buffer[1]='\0';
	  return buffer;
	}
	
	for (i=0; value > 0; i++, value /= radix)
	  buffer[i] = base[(value % radix)];
	
	buffer[i] ='\0';

	buffer = strreverse(buffer, i);
	return buffer;
}



/*-------------------------------- attachext ---------------------------------
 *    
 * attaches an extension to a filename
 * 
 */
 
char *
attachext (void *space, char *str, Uint l, char *ext, Uint m)
{
    char *new;

	new = ALLOCMEMORY(space, NULL, char, l+m+1);
	strncpy(new, str, l);
	new[l]='\0';
	strncat(new, ext, m);
	
	return new;
}
/*-------------------------------- attachpath --------------------------------
 *    
 * attach a path (or any other string) to a filename (or any other string) 
 * and extensions (or any other string) at the end.
 * 
 */
 
char*
attachpath (void *space, char *str, Uint l, char *path, Uint m,  
				char *ext, Uint n)
{
	char *new;
	
	new = ALLOCMEMORY(space, NULL, char, (l+m+n+1));
	strncpy(new, path, m);
	new[m] = '\0';
	strncat(new, str, l);
	strncat(new, ext, n);
	
	return new;
}

int checkmd5(unsigned char *a, unsigned char *b) { 
  return (strncmp((char*)a, (char*)b, 16));
}

void
fprintStringset(FILE *dev, stringset_t *set) {
  Uint i;

  for(i=0; i < set->noofstrings; i++) {
   fprintf(dev, "%d:'%s' (len:%d)\n", i, 
       set->strings[i].str, set->strings[i].len);
  }

}


char *
sprintchar (char **str, char chr) {
  char *tmp;

  tmp = calloc(2, 1);
  tmp[0] = chr;

  *str = tmp;
  return *str;
}

char *
sprintint (char **str, int n) {
  char *tmp;

  tmp = calloc(MAX_INT_LENGTH+1, 1);
  sprintf(tmp, "%d", n);

  *str = tmp;
  return *str;
}

char *
sprintUint (char **str, Uint n) {
  char *tmp;

  tmp = calloc(MAX_INT_LENGTH+1, 1);
  sprintf(tmp, "%u", n);

  *str = tmp;
  return *str;
}

char *
sprintstr (char **str, char *src, Uint len) {
  char *tmp;

  tmp = calloc(len+1, 1);
  memmove(tmp, src, len);

  *str = tmp;
  return *str;
}

char *
sprintflt (char **str, double flt) {
  Uint size;
  char *tmp;

  size = snprintf(NULL, 0, "%.4f", flt);
  tmp = calloc(size+2, 1);
  snprintf(tmp, size+1, "%.4f", flt);

  *str = tmp;
  return *str;
}

char *my_strdup(const char *str) {
  size_t len = strlen(str);
  char *x = (char *)malloc(len+1); /* 1 for the null terminator */
  if(!x) return NULL; /* malloc could not allocate memory */
  memcpy(x,str,len); /* copy the string into the new buffer */
  x[len]=0;
  return x;
}
