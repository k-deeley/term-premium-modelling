function mustBeSquare( X )
%MUSTBESQUARE Validate that the input, X, is square.

assert( ismatrix( X ), "mustBeSquare:Not2D", ...
    "The input is not 2-D." )
assert( width( X ) == height( X ), "mustBeSquare:NotSquare", ...
    "The input is not square." )

end % mustBeSquare