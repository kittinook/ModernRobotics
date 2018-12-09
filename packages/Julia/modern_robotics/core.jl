from __future__ import print_function
'''
***************************************************************************
Modern Robotics: Mechanics, Planning, and Control.
Code Library
***************************************************************************
Author: Mark Saroufim
Email: marksaroufim@gmail.com
Date: December 2018
***************************************************************************
Language: Julia
Also available in: Python, MATLAB, Mathematica
Required library: LinearAlgebra
Optional library:
***************************************************************************
'''

'''
*** IMPORTS ***
'''
using LinearAlgebra
'''
*** BASIC HELPER FUNCTIONS ***
'''
Vector{T} = Array{T,1}

function NearZero(z::AbstractFloat)::Bool
    """Determines whether a scalar is small enough to be treated as zero

    :param z: A scalar input to check
    :return: True if z is close to zero, false otherwise

    Example Input:
        z = -1e-7
    Output:
        true
    """
    return abs(z) < 1e-6
end

function Normalize(V::Vector)::Vector
    """Normalizes a vector

    :param V: A vector
    :return: A unit vector pointing in the same direction as z

    Example Input:
        V = [1 2 3])
    Output:
        np.array([0.26726124 0.53452248 0.80178373])
    """
    return normalize(V)
end

'''
*** CHAPTER 3: RIGID-BODY MOTIONS ***
'''

function RotInv(R::Matrix)::Matrix
    """Inverts a rotation matrix

    :param R: A rotation matrix
    :return: The inverse of R

    Example Input:
        R = [[0, 0, 1],
             [1, 0, 0],
             [0, 1, 0]]
    Output:
        ([[0, 1, 0],
          [0, 0, 1],
          [1, 0, 0]])
    """
    return transpose(R)
end

#should be a cleaner way to subclass this from Vector and still have named arguments
struct Vec3
    x::AbstractFloat
    y::AbstractFloat
    z::AbstractFloat
end

mutable struct so3
    S::Array{Float64,2}
    so3() = new()
    so3(s::FloatInt) = new(zeros(3,3))
    so3(v::VectorFloatInt) = new(skew(v))
    so3(S::Array{Float64,2}) = new(S)
end

function VecToso3(omg::Vec3)::so3
    """Converts a 3-vector to an so(3) representation

    :param omg: A 3-vector
    :return: The skew symmetric representation of omg

    Example Input:
        omg = [1, 2, 3]
    Output:
        [ 0, -3,  2],
        [ 3,  0, -1],
        [-2,  1,  0]])
    """
    return [0      -omg[3]  omg[2];
            omg[3]       0 -omg[1];
            -omg[2] omg[2]  0;]
end

function so3ToVec(so3mat::so3):
    """Converts an so(3) representation to a 3-vector

    :param so3mat: A 3x3 skew-symmetric matrix
    :return: The 3-vector corresponding to so3mat

    Example Input:
        so3mat = [[ 0, -3,  2],
                  [ 3,  0, -1],
                  [-2,  1,  0]])
    Output:
        np.array([1, 2, 3])
    """
    return [so3mat[3][2] so3mat[1][3] so3mat[2][1]]
end

# think about a good return type for this
function AxisAng3(expc3)
    return (normalize(expc3), norm(expc3))
end

function eye(x)::Matrix
    return Diagonal(ones(x,x))
end


function MatrixExp3(so3mat::so3)
    omgtheta = so3toVec(so3mat)
    if NearZero(norm(omgtheta))
        return eye(3,3)
    else
        theta = AxisAng3(omgtheta)[1]
        omgmat = so3mat / theta
        return eye(3,3)
    end
end
