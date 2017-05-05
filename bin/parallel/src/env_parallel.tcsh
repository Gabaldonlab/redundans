#!/usr/bin/env tcsh

# This file must be sourced in tcsh:
#
#   source `which env_parallel.tcsh`
#
# after which 'env_parallel' works
#
#
# Copyright (C) 2016,2017
# Ole Tange and Free Software Foundation, Inc.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
# or write to the Free Software Foundation, Inc., 51 Franklin St,
# Fifth Floor, Boston, MA 02110-1301 USA

set _parallel_exit_CODE=0
if ("`alias env_parallel`" == '') then
  # Activate alias
  alias env_parallel '(setenv PARALLEL "\!*"; source `which env_parallel.tcsh`)'
else
  # Get the --env variables if set
  # --env _ should be ignored
  # and convert  a b c  to (a|b|c)
  # If --env not set: Match everything (.*)
  set _tMpscRIpt=`tempfile`
  cat <<'EOF' > $_tMpscRIpt
            #!/usr/bin/perl

            for(@ARGV){
                /^_$/ and $next_is_env = 0;
                $next_is_env and push @envvar, split/,/, $_;
                $next_is_env = /^--env$/;
            }
            $vars = join "|",map { quotemeta $_ } @envvar;
            print $vars ? "($vars)" : "(.*)";
'EOF'
  set _grep_REGEXP="`perl $_tMpscRIpt -- $PARALLEL`"

  # Deal with --env _
  cat <<'EOF' > $_tMpscRIpt
            #!/usr/bin/perl
            
            for(@ARGV){
                $next_is_env and push @envvar, split/,/, $_;
                $next_is_env=/^--env$/;
            }
            if(grep { /^_$/ } @envvar) {
                if(not open(IN, "<", "$ENV{HOME}/.parallel/ignored_vars")) {
             	    print STDERR "parallel: Error: ",
            	    "Run \"parallel --record-env\" in a clean environment first.\n";
                } else {
            	    chomp(@ignored_vars = <IN>);
            	    $vars = join "|",map { quotemeta $_ } "env_parallel", @ignored_vars;
		    print $vars ? "($vars)" : "(,,nO,,VaRs,,)";
                }
            }
'EOF'
  set _ignore_UNDERSCORE="`perl $_tMpscRIpt -- $PARALLEL`"
  rm $_tMpscRIpt

  # Get the scalar and array variable names
  set _vARnAmES=(`set | awk -e '{print $1}' |grep -vE '^(#|_|killring|prompt2|command)$' | grep -Ev '^(PARALLEL_TMP)$' | grep -E "^$_grep_REGEXP"\$ | grep -vE "^$_ignore_UNDERSCORE"\$`)

  # Make a tmpfile for the variable definitions
  set _tMpvARfILe=`tempfile`
  
  # Make a tmpfile for the variable definitions + alias
  set _tMpaLLfILe=`tempfile`
  foreach _vARnAmE ($_vARnAmES);
    # if not defined: next
    eval if'(! $?'$_vARnAmE') continue'
    # if $#myvar <= 1 echo scalar_myvar=$var
    eval if'(${#'$_vARnAmE'} <= 1) echo scalar_'$_vARnAmE'='\"\$$_vARnAmE\" >> $_tMpvARfILe;
    # if $#myvar > 1 echo array_myvar=$var
    eval if'(${#'$_vARnAmE'} > 1) echo array_'$_vARnAmE'="$'$_vARnAmE'"' >> $_tMpvARfILe;
  end
  unset _vARnAmE _vARnAmES
  # shell quote variables (--plain needed due to $PARALLEL abuse)
  # Convert 'scalar_myvar=...' to 'set myvar=...'
  # Convert 'array_myvar=...' to 'set array=(...)'
  cat $_tMpvARfILe | parallel --plain --shellquote |  perl -pe 's/^scalar_(\S+).=/set $1=/ or s/^array_(\S+).=(.*)/set $1=($2)/ && s/\\ / /g;' > $_tMpaLLfILe
  # Cleanup
  rm $_tMpvARfILe; unset _tMpvARfILe

# ALIAS TO EXPORT ALIASES:

#   Quote ' by putting it inside "
#   s/'/'"'"'/g;
#   ' => \047 " => \042
#   s/\047/\047\042\047\042\047/g;
#   Quoted: s/\\047/\\047\\042\\047\\042\\047/g\;

#   Remove () from second column
#   s/^(\S+)(\s+)\((.*)\)/\1\2\3/;
#   Quoted: s/\^\(\\S+\)\(\\s+\)\\\(\(.\*\)\\\)/\\1\\2\\3/\;

#   Add ' around second column
#   s/^(\S+)(\s+)(.*)/\1\2'\3'/
#   \047 => '
#   s/^(\S+)(\s+)(.*)/\1\2\047\3\047/;
#   Quoted: s/\^\(\\S+\)\(\\s+\)\(.\*\)/\\1\\2\\047\\3\\047/\;

#   Quote ! as \!
#   s/\!/\\\!/g;
#   Quoted: s/\\\!/\\\\\\\!/g;

#   Prepend with "\nalias "
#   s/^/\001alias /;
#   Quoted: s/\^/\\001alias\ /\;
  alias | \
    grep -E "^$_grep_REGEXP" | \
    grep -vE "^$_ignore_UNDERSCORE""[^_a-zA-Z]" | \
    perl -pe s/\\047/\\047\\042\\047\\042\\047/g\;s/\^\(\\S+\)\(\\s+\)\\\(\(.\*\)\\\)/\\1\\2\\3/\;s/\^\(\\S+\)\(\\s+\)\(.\*\)/\\1\\2\\047\\3\\047/\;s/\^/\\001alias\ /\;s/\\\!/\\\\\\\!/g >> $_tMpaLLfILe 
  
  setenv PARALLEL_ENV "`cat $_tMpaLLfILe; rm $_tMpaLLfILe`";
  unset _tMpaLLfILe;
  # Use $PARALLEL set in calling alias
  parallel
  set _parallel_exit_CODE=$status
  setenv PARALLEL_ENV
  setenv PARALLEL
endif
(exit $_parallel_exit_CODE)

# Tested working for aliases
# alias env_parallel 'setenv PARALLEL_ENV "`alias | perl -pe s/\\047/\\047\\042\\047\\042\\047/g\;s/\^\(\\S+\)\(\\s+\)\\\(\(.\*\)\\\)/\\1\\2\\3/\;s/\^\(\\S+\)\(\\s+\)\(.\*\)/\\1\\2\\047\\3\\047/\;s/\^/\\001alias\ /\;s/\\\!/\\\\\\\!/g;`";parallel \!*; setenv PARALLEL_ENV'

