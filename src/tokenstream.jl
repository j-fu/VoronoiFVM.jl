mutable struct TokenStream
    input::IOStream
    tokens::Array{SubString{String},1}
    itoken::Int
    lineno::Int
    comment::Char
    dlm::Function
    TokenStream()=new()
end

destruct!(tks::TokenStream)=close(tks.input)

struct UnexpectedTokenError <: Exception
    found::String
    expected::String
    lineno::Int
end

Base.showerror(io::IO, e::UnexpectedTokenError) = print(io,"Unexpected token in line $(e.lineno): $(e.found) (expected $(e.expected))")

export TokenStream,gettoken, expecttoken,trytoken

"""
Create Tokenstream.
"""
TokenStream(filename::String;comment='#', dlm=isspace)=TokenStream(open(filename),comment=comment,dlm=dlm)

"""
Create Tokenstream.
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
Get next token.
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
  Try for keyword token.
  It token is missing, the token read is put back into stream
  and the next try/gettoken command continues at the same position.
  A value of false is returned.
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

