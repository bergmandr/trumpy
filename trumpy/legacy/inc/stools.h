/*
 *  Tools for manipulating strings that may or may not already exist.
 */

/*
 *  returns '0' if 'str' does not end with 'sfx', '1' if it does, and '-1' if 
 *    'sfx' is longer than 'str'.
 */
int endswith(char *str, const char *sfx);

/*
 *  returns '0' if 'str' does not start with 'pfx', '1' if it does, and '-1' 
 *    if 'pfx' is longer than 'str'.
 */
int startswith(char *str, const char *pfx);
