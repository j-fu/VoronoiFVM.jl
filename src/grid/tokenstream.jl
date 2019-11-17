#=
Token stream module
=#

"""
$(TYPEDEF)

Tokenstream allows to read tokenized data from file without keeping
the file ocntent in memory.

$(TYPEDFIELDS)
"""
mutable struct TokenStream
    """
    Input stream
    """
    input::IOStream

    """
    Array of current tokens kept in memory.
    """
    tokens::Array{SubString{String},1}

    """
    Position of actual token in tokens array
    """
    itoken::Int


    """
    Line number in IOStream
    """
    lineno::Int

    """
    Comment character
    """
    comment::Char

    """
    Function telling if given character is a delimiter.
    """
    dlm::Function
    TokenStream()=new()
end

"""
  $(TYPEDSIGNATURES) 
    
    Tokenstream destructor should close input
"""
destruct!(tks::TokenStream)=close(tks.input)


"""
$(TYPEDEF)

Error thrown when the token expected  in expect!  is not there.
$(TYPEDFIELDS)
"""
struct UnexpectedTokenError <: Exception
    found::String
    expected::String
    lineno::Int
end

Base.showerror(io::IO, e::UnexpectedTokenError) = print(io,"Unexpected token in line $(e.lineno): $(e.found) (expected $(e.expected))")


"""
  $(TYPEDSIGNATURES) 
        
Create Tokenstream with file name argument.
"""
TokenStream(filename::String;comment='#', dlm=isspace)=TokenStream(open(filename),comment=comment,dlm=dlm)

"""
  $(TYPEDSIGNATURES) 

Create Tokenstream with IOStream argument.
"""
function TokenStream(input::IOStream; comment='#', dlm=isspace)
    tks=TokenStream()
    tks.input=input
    tks.lineno=0
    tks.comment=comment
    tks.dlm=dlm
    tks.tokens=split("")
    _fetch!(tks)
    return tks
end

"""
  $(TYPEDSIGNATURES) 

Check if all tokens have been consumed.
"""
Base.eof(tks::TokenStream)=!_fetch!(tks)

#=
Read next line with tokens and split it.
Return if tokens are left
=#
function _fetch!(tks::TokenStream)
    if (tks.itoken<=length(tks.tokens))
        return true
    end
    hasline=false
    buffer=""
    while !hasline && !eof(tks.input)
        buffer=readline(tks.input)
        tks.lineno+=1
        if (length(buffer)>0) && (buffer[1]!=tks.comment)
            hasline=true
        end
    end
    tks.tokens=split(buffer,tks.dlm,keepempty=false)
    tks.itoken=1
    return (tks.itoken<=length(tks.tokens))
end


"""
  $(TYPEDSIGNATURES) 

Get next token from tokenstream.
"""
function  gettoken(tks::TokenStream)
    if _fetch!(tks)
        token=tks.tokens[tks.itoken]
        tks.itoken+=1
        return token
    end
    nothing
end

"""
  $(TYPEDSIGNATURES) 

  Expect keyword token.

  If token is missing, an UnexpectedTokenError is thrown
  If the token  has been found, reading will continue  at the position
  after the token found.
"""
function expecttoken(tks::TokenStream,expected::String)
    token=gettoken(tks)
    if token!=expected
        throw(UnexpectedTokenError(token,expected,tks.lineno))
    end
    true
end

"""
  $(TYPEDSIGNATURES) 

  Try for keyword token.


  It token is missing, the token read is put back into stream,
  a value of false is returned and the next try/gettoken command continues at the same position,

  Otherwise, true is returned, and reading continues after the token found.
"""
function trytoken(tks::TokenStream,expected::String)
    token=gettoken(tks)
    if token!=expected
        tks.itoken=tks.itoken-1
        return false
    end
    return true
end

