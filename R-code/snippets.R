
# convert all columns of one type to another
DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], as.factor)


# LaTeX: toggle whether or not to display a section

\newif\ifdraft
\drafttrue % or \draftfalse

\ifdraft
% text
\ifelse
%alternative text
\fi