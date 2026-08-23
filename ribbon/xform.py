#!/usr/bin/env python
"""
**Euler angle ranges**

=============   ===============   ================  ================
 Sequence            phi               theta           psi
=============   ===============   ================  ================
  XYZ,   ZYX     [-pi,   pi]        [-pi/2, pi/2]     [-pi, pi]
  XZY,   YZX     [-pi,   pi]        [-pi,  pi]        [-pi/2, pi/2]
  ZXY,   YXZ     [-pi/2, pi/2]      [-pi,  pi]        [-pi, pi]
=============   ===============   ================  ================

Euler angle sequence: 'XYZ' (world). First rotation by *phi* about X, second
rotation by *theta* about Y, and the third rotation by *psi* about Z axis of
the world (i.e. fixed) frame.

This is the same as the sequence used in the `Blender blenlib code
<https://github.com/blender/blender/blob/
7d641fe96810cdc2598b2f37ec4f6192e133e6d8/
source/blender/blenlib/BLI_math_euler_types.hh>`__.

In contrast, the 'XYZ' sequence is understood in the Aerospace community as:
First rotation about Z-axis, second rotation about Y-axis, and the third
rotation about X-axis of the body frame.

*Reference:* D. Eberly, `Euler angle formulas.
<http://www.geometrictools.com/Documentation/EulerAngles.pdf>`__

"""

import math
import numpy as np

def translate(v, delta):
    '''
    Translates vectors in-place by delta.

    Parameters
    ----------
    v : (n,3) ndarray
        Vectors to translate.
    delta : (3,) ndarray
        Translation vector.

    Returns
    -------
    (n,3) ndarray
        Translated vectors
    '''
    n = v.shape[0]
    for i in range(n):
        v[i,:] += delta
    return v


def align(v, old, new):
    '''
    Rotate vectors `v` such that `old` and `new` align.
    represent coordinate axes. They must be unit vectors.

    Parameters
    ----------
    v : (n,3) ndarray
        Vectors to rotate.
    old : ndarray
        Coordinate axes.
    new : ndarray
        Coordinate axes.

    Returns
    -------
    (n,3) ndarray
        Rotated vectors.
    '''
    assert old.shape[0] == new.shape[0]
    n = old.ndim
    if n == 1:
        #Angle between old and new
        cos_angle = np.dot(old, new)
        if (cos_angle > 1.0) or (cos_angle < -1.0):
            if math.isclose(abs(cos_angle), 1.0):
                cos_angle = math.copysign(1.0, cos_angle)
            else:
                mag_old = np.linalg.norm(old)
                mag_new = np.linalg.norm(new)
                raise ValueError(
                    f"Both `old`(= {old}) and `new`(= {new}) must be unit"
                    f" vectors. Magnitude of `old` = {mag_old} and"
                    f" `new` = {mag_new}.")
        angle = math.acos(cos_angle)
        #Axis of rotation 
        axis = np.cross(old, new)
        axis_nrm = np.linalg.norm(axis)
        if math.isclose(axis_nrm, 0.0, rel_tol=1e-14):
            raise ValueError(f"`old`(= {old}) and `new`(= {new}) are collinear.")
        else:
            axis /= axis_nrm
        return rotate_vector_axis_angle(v, axis, angle)

    elif n == 2:
        z_old = np.cross(old[0,:], old[1,:])
        z_old_nrm = np.linalg.norm(z_old)
        if math.isclose(z_old_nrm, 0.0, rel_tol=1e-14):
            raise ValueError( f"Axes 0 and 1 of `old` are collinear.\n"
                             f"Axis 0 = {old[0,:]}.\nAxis 1 = {old[1,:]}")
        else:
            z_old /= z_old_nrm

        z_new = np.cross(new[0,:], new[1,:])
        z_new_nrm = np.linalg.norm(z_new)
        if math.isclose(z_new_nrm, 0.0, rel_tol=1e-14):
            raise ValueError( f"Axes 0 and 1 of `new` are collinear.\n"
                             f"Axis 0 = {new[0,:]}.\nAxis 1 = {new[1,:]}")
        else:
            z_new /= z_new_nrm

        axes_old = np.vstack((old, z_old))
        axes_new = np.vstack((new, z_new))
        dcm = dcm_from_axes(axes_old, axes_new)
        return rotate_vector_dcm(v, dcm)
    elif n == 3:
        dcm = dcm_from_axes(old, new)
        return rotate_vector_dcm(v, dcm)


def mat_is_rotmat(mat):
    """
    Checks if `mat` is a rotation matrix or not.

    Parameters
    ----------
    mat : (3,3) ndarray
        Array to check.

    Returns
    -------
    bool
        ``True`` if `mat` is a rotation matrix, ``False`` otherwise.
    """
    det_is_one = math.isclose(np.linalg.det(mat), 1.0, 
                              abs_tol=1e-12, rel_tol=1e-12)
    is_orthogonal = np.allclose(np.dot(mat, mat.T), np.identity(3))
    return is_orthogonal and det_is_one


def get_ransphere(rad, rng=None):
    """
    Generates a random vector distributed uniformly on the surface of a sphere.

    Parameters
    ----------
    rad : float
        Radius of the sphere.
    rng : :py:class:`numpy.random.Generator`
        A random number generator. If ``None``, the generator is constructed
        with :py:func:`numpy.random.default_rng()`.

    Returns
    -------
    (3,) ndarray
        Random vector of magnitude `rad`.

    """
    rng_ = rng if rng else np.random.default_rng()
    ransphere = np.zeros((3,))
    while True:
        zeta1 = -1.0 + 2*rng_.random()
        zeta2 = -1.0 + 2*rng_.random()
        zetasq = zeta1*zeta1 + zeta2*zeta2
        if zetasq <= 1.0:
            rt = math.sqrt(1.0 - zetasq)
            ransphere[0] = 2*zeta1*rt
            ransphere[1] = 2*zeta2*rt
            ransphere[2] = 1.0 - 2*zetasq
            break
    ransphere = ransphere/np.linalg.norm(ransphere)
    return rad*ransphere


def get_frc_link(ri, theta, rng=None):
    """
    Given a segment `ri` of a freely rotating chain, return the unit vector `r`
    along the next link such that *r(i+1) = ri + l \* r*, where *l* is the length of
    the link.

    Parameters
    ----------
    ri : (3,) ndarray
        Vector representing link `i`. 
    theta : float
        Rotation angle between the links in radian.
    rng : :py:class:`numpy.random.Generator`
        A random number generator. If ``None``, the generator is constructed
        with :py:func:`numpy.random.default_rng()`.

    Returns
    -------
    (3,) ndarray
        Unit vector along link *i+1*.
    """
    rng_ = rng if rng else np.random.default_rng()
    phi = 2*math.pi*rng_.random()
    ctheta = math.cos(theta)
    stheta = math.sin(theta)
    cphi = math.cos(phi)
    sphi = math.sin(phi)
    
    rihat = ri/np.linalg.norm(ri)
    
    yhat = np.array([0, 1, 0])
    axis = np.cross(yhat, rihat)
    axis = axis/np.linalg.norm(axis)
    cos_angle = np.dot(rihat, yhat)
    angle = math.acos(cos_angle)
    axis, angle = fix_axis_angle(axis, angle)
    
    r = np.array([stheta*cphi, ctheta, stheta*sphi])
    r = rotate_vector_axis_angle(r, axis, angle)
    return r


#QUATERNION-----------------------------------------------------------

def get_rand_quat(rng=None):
    """
    Returns a random unit quaternion.

    Parameters
    ----------
    rng : :py:class:`numpy.random.Generator`
        A random number generator. If ``None``, the generator is constructed
        with :py:func:`numpy.random.default_rng()`.

    Returns
    -------
    (4,) ndarray
        Unit quaternion.
    """
    rng_ = rng if rng else np.random.default_rng()
    q = rng_.random((4,))
    return normalize_quat(q)


def get_identity_quat():
    """
    Returns the identity unit quaternion.

    Returns
    -------
    (4,) ndarray
        Unit quaternion.
    """
    return np.array([1.0, 0.0, 0.0, 0.0])


def conjugate_quat(q):
    '''
    Conjugates a quaternion in-place and returns it.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Conjugated quaternion.
    '''
    q[1:4] = -q[1:4]
    return q


def get_conjugated_quat(q):
    """
    Returns a conjugated quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Conjugated quaternion.
    """
    p = np.copy(q)
    p[1:4] = -p[1:4]
    return p


def invert_quat(q):
    '''
    Inverts a quaternion in-place and returns it.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Inverted quaternion.
    '''
    return conjugate_quat(q)


def get_inverted_quat(q):
    '''
    Returns an inverted quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Inverted quaternion.    
    '''
    p = np.copy(q)
    return conjugate_quat(p)


def normalize_quat(q):
    '''
    Normalizes a quaternion in-place and returns it.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Normalized quaternion.
    '''
    q /= np.linalg.norm(q)
    return q


def get_normalized_quat(q):
    '''
    Returns a normalized quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    (4,) ndarray
        Normalized quaternion.
    '''
    p = np.copy(q)
    return normalize_quat(p)


def quat_is_normalized(q):
    """
    Checks whether a quaternion is normalized or not.

    Parameters
    ----------
    q : (4,) ndarray
        Quaternion.

    Returns
    -------
    bool
        ``True`` if `q` is normalized, ``False`` otherwise.
    """
    norm = np.linalg.norm(q)
    return math.isclose(norm, 1.0, rel_tol=1e-14)


def get_quat_prod(p, q):
    """
    Returns the product of two unit quaternions.

    Parameters
    ----------
    p : (4,) ndarray
        Unit quaternion.
    q : (4,) ndarray
        Unit quaternion.
    
    Returns
    -------
    (4,) ndarray
        Product of two unit quaternions. This is also a unit quaternion.
    
    """
    p0, p1, p2, p3 = tuple(p)
    prod_mat = np.array([[p0, -p1, -p2, -p3],
                         [p1,  p0, -p3,  p2],
                         [p2,  p3,  p0, -p1],
                         [p3, -p2,  p1,  p0]])
    pq = normalize_quat(np.dot(prod_mat, q))
    return pq


def interpolate_quat(q1, q2, t):
    """
    Interpolate between two unit quaternions.

    Parameters
    ----------
    q1 : (4,) ndarray
        Unit quaternion.
    q2 : (4,) ndarray
        Unit quaternion.
    t : float
        A fraction between 0 and 1 (both inclusive) specifying the
        interpolation point.
    
    Returns
    -------
    (4,) ndarray
        Interpolated unit quaternion.
    """
    theta = get_angle_between_quat(q1, q2)
    q = (q1*math.sin((1.0-t)*theta)
            + q2*math.sin(t*theta))/math.sin(theta)
    return normalize_quat(q)


def get_angle_between_quat(p, q):
    """
    Returns the angle between two unit quaternions p and q.

    Parameters
    ----------
    p : (4,) ndarray
        Unit quaternion.
    q : (4,) ndarray
        Unit quaternion.
    
    Returns
    -------
    float
        Angle in radian.
    """
    return math.acos(np.dot(p,q))


def quat_deriv_to_ang_vel(q, qdot):
    """
    Calculates the angle velocity from a unit quaternion and its time
    derivative.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    qdot : (4,) ndarray
        Derivative of `q`.
    
    Returns
    -------
    (3,) ndarray
        Angular velocity.
    """
    mat = quat_deriv_to_ang_vel_mat(q)
    return np.dot(mat, qdot)


def quat_deriv_to_ang_vel_mat(qdot):
    """
    Returns the matrix mapping the derivative of a unit quaternion to angular
    velocity.

    Parameters
    ----------
    qdot : (4,) ndarray
        Derivative of a unit quaternion.
    
    Returns
    -------
    (3,4) ndarray
        Angular velocity matrix.
    """
    q0, q1, q2, q3 = tuple(qdot)
    return 2*np.array([[-q1,  q0, -q3,  q2],
                       [-q2,  q3,  q0, -q1],
                       [-q3, -q2,  q1,  q0]])


def ang_vel_to_quat_deriv(q, ang_vel):
    """
    Calculates the time derivative of a unit quaternion from angular velocity.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    ang_vel : (3,) ndarray
        Angular velocity.
    
    Returns
    -------
    (3,) ndarray
        Derivative of `q`.
    """
    mat = ang_vel_to_quat_deriv_mat(q)
    qdot = np.dot(mat, ang_vel)
    return qdot


def ang_vel_to_quat_deriv_mat(q):
    """
    Returns the matrix mapping angular velocity to time derivative of a unit
    quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    (4,3) ndarray
        Quaternion derivative matrix.
    """
    q0, q1, q2, q3 = tuple(q)
    return 0.5*np.array([[-q1, -q2, -q3],
                         [ q0,  q3, -q2],
                         [-q3,  q0,  q1],
                         [ q2, -q1,  q0]])


def quat_to_axis_angle(q):
    """
    Converts a unit quaternion to an *axis-angle* representation.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    axis : (3,) ndarray
        Unit vector along the axis of rotation.
    angle : float
        Angle in radian.
    """
    angle = 2*math.acos(q[0])
    sin = math.sqrt(1.0-q[0]**2)
    if angle > 0.0:
        if angle < math.pi:
            axis = q[1:4]/sin
        else:
            rotmat = get_rotmat_quat(q)
            axis, angle = extract_axis_angle_from_rotmat(rotmat)
    else:
        axis = np.array([1.0, 0.0, 0.0])
    return fix_axis_angle(axis, angle, normalize=True)


def quat_to_euler(q, seq='XYZ', world=True):
    """
    Converts a unit quaternion to an Euler angle sequence.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.
        seq

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = get_rotmat_quat(q)
    return factorize_rotmat(rotmat, seq=seq, world=world)


def quat_to_dcm(q):
    """
    Converts a unit quaternion to a direction cosine matrix.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix.
    """
    return get_shiftmat_quat(q, forward=True)


def any_to_quat(orientation):
    """Converts an orientation to unit quaternion from either of the
    following: (1) Quaternion, (2) Euler angles, (3) Axis-angle, or (4)
    Direction cosine matrix.

    Parameters
    ----------
    orientation : dict
        The keys and values are:

        - `'repr'`  = ``'quat'`` | ``'euler'`` | ``'axis_angle'`` | ``'dcm'``
        - `'quat'`  = (4,) *ndarray*
        - `'euler'` = (3,) *ndarray*
        - `'seq'` = ``'XYZ'`` | ``'XZY'`` | ``'YXZ'`` | ``'YZX'``
          | ``'ZXY'`` | ``'ZYX'``
        - `'world'` = ``True`` | ``False``
        - `'axis'` = (3,) *ndarray*
        - `'angle'` = *float*
        - `'dcm'` = (3,3) *ndarray*

        For any value of `'repr'`, only the relevant keys are accessed, the
        rest are ignored. E.g., if `'repr'` = ``'quat'``, only the `'quat'` key
        is necessary, but for `'repr'` = ``'euler'``, the required keys are
        `'euler'`, `'seq'`, and `'world'`.

    Returns
    -------
    (4,) ndarray
        Unit quaternion.
    """
    ori_repr = orientation['repr']
    if ori_repr == 'quat':
        quat = np.array(orientation['quat'])
    elif ori_repr == 'euler':
        euler = np.array(orientation['euler'])
        seq = orientation['seq']
        world = orientation['world']
        quat = euler_to_quat(euler, seq=seq, world=world)
    elif ori_repr == 'axis_angle':
        axis = np.array(orientation['axis'])
        angle = orientation['angle']
        quat = axis_angle_to_quat(axis, angle)
    elif ori_repr == 'dcm':
        quat = dcm_to_quat(orientation['dcm'])
    else:
        raise ValueError(
            'Unrecognized orientation repr {0}'.format(ori_repr))
    return quat


def rotate_vector_quat(v, q):
    """
    Rotates vectors by a unit quaternion.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to rotate.
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.
    """
    rotmat = get_rotmat_quat(q)
    return np.dot(v, rotmat.T)


def get_rotmat_quat(q):
    """Returns the rotation matrix corresponding to a unit quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix.
    """
    rotmat = np.empty((3,3))

    q0sq = q[0]*q[0]
    q1sq = q[1]*q[1]
    q2sq = q[2]*q[2]
    q3sq = q[3]*q[3]
    q0q1 = q[0]*q[1]
    q0q2 = q[0]*q[2]
    q0q3 = q[0]*q[3]
    q1q2 = q[1]*q[2]
    q1q3 = q[1]*q[3]
    q2q3 = q[2]*q[3]

    rotmat[0,0] = 2*(q0sq + q1sq) - 1.0
    rotmat[0,1] = 2*(q1q2 - q0q3)
    rotmat[0,2] = 2*(q1q3 + q0q2)
    rotmat[1,0] = 2*(q1q2 + q0q3)
    rotmat[1,1] = 2*(q0sq + q2sq) - 1.0
    rotmat[1,2] = 2*(q2q3 - q0q1)
    rotmat[2,0] = 2*(q1q3 - q0q2)
    rotmat[2,1] = 2*(q2q3 + q0q1)
    rotmat[2,2] = 2*(q0sq + q3sq) - 1.0
    return rotmat


def shift_vector_quat(v, q, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a unit quaternion, shifts vectors from A
    to B or B to A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to shift.
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,) or (n,3) ndarray
        Shifted vectors.
    """
    shiftmat = get_shiftmat_quat(q, forward=forward)
    return np.dot(v, shiftmat.T)


def shift_tensor2_quat(a, quat, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a unit quaternion, shifts second-order tensors from A
    to B or B to A.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor.
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifted second-order tensor.
    """
    shiftmat = get_shiftmat_quat(quat, forward=forward)
    return np.einsum('ip,jq,pq', shiftmat, shiftmat, a)


def shift_tensor3_quat(a, quat, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a unit quaternion, shifts third-order tensors from A
    to B or B to A.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor.
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3,3) ndarray
        Shifted third-order tensor.
    """
    shiftmat = get_shiftmat_quat(quat, forward=forward)
    return np.einsum('ip,jq,kr,pqr', shiftmat, shiftmat, shiftmat, a)


def get_shiftmat_quat(q, forward=False):
    """
    Returns the shifter matrix corresponding to a unit quaternion.

    Parameters
    ----------
    q : (4,) ndarray
        Unit quaternion.
    forward : bool
       Whether to shift forward, i.e., along the orientation or shift reverse.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.
    """
    if forward:
        shiftmat = get_rotmat_quat(get_conjugated_quat(q))
    else:
        shiftmat = get_rotmat_quat(q)
    return shiftmat


#AXIS-ANGLE------------------------------------------------------------
def fix_axis_angle(axis, angle, normalize=True):
    """
    Returns a copy of `axis` and `angle` by modifying their values such that it
    is a right handed rotation with `angle` lies in [0, *pi*].

    Parameters
    ----------
    axis : (3,) ndarray
        Axis of rotation (need not be a unit vector).
    angle : float
        Angle in radian.
    normalize : bool
        Whether to normalize the axis to a unit vector.

    Returns
    -------
    axis : (3,) ndarray
        Modified axis of rotation, possibly normalized.
    angle : float
        Modified angle in radian.

    """
    if normalize:
        norm = np.linalg.norm(axis)
        if not math.isclose(norm, 1.0, abs_tol=1e-14, rel_tol=1e-14):
            axis /= norm
    angle = math.fmod(angle, 2*math.pi)
    if angle < 0.0:
        angle = -angle
        axis = -axis
    if angle > math.pi:
        angle = 2*math.pi - angle
        axis = -axis
    return (axis, angle)


def get_rand_axis_angle():
    """
    Generates a random orientation in *axis-angle* representation.

    The axis is a random vector drawn from a uniform distribution on the
    surface of a unit sphere. The current implementation in based on the
    algorithm from Allen & Tildesley p. 349.

    Returns
    -------
    axis : (3,) ndarray
        Axis of rotation. This is a unit vector.
    angle : float
        Angle in radian.
    """
    axis = np.zeros((3,))
    #Generate angle: A uniform random number from [0.0, 2*pi)
    angle = 2.0*math.pi*np.random.random()
    while True:
        #Generate two uniform random numbers from [-1, 1)
        zeta1 = 2.0*np.random.random() - 1.0
        zeta2 = 2.0*np.random.random() - 1.0
        zetasq = zeta1**2 + zeta2**2
        if zetasq <= 1.0:
            break
    rt = np.sqrt(1.0-zetasq)
    axis[0] = 2.0*zeta1*rt
    axis[1] = 2.0*zeta2*rt
    axis[2] = 1.0 - 2.0*zetasq
    return fix_axis_angle(axis, angle)


def axis_angle_to_quat(axis, angle):
    """
    Converts an *axis-angle* representation to a unit quaternion.

    Parameters
    ----------
    axis : (3,) ndarray
        Axis of rotation. This must be a unit vector.
    angle : float
        Angle

    Returns
    -------
    q : (4,) ndarray
        Unit quaternion.

    """
    w = math.cos(angle/2)
    v = math.sin(angle/2)*axis
    q = np.array([w, v[0], v[1], v[2]])
    return normalize_quat(q)


def axis_angle_to_euler(axis, angle, seq='XYZ', world=True):
    """Coverts an *axis-angle* representation to *Euler angles*.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the direction of the axis.
    angle : float
        Angle in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = get_rotmat_axis_angle(axis, angle)
    euler = factorize_rotmat(rotmat, seq=seq, world=world)
    return euler


def axis_angle_to_dcm(axis, angle):
    """Converts an *axis-angle* representation to a direction cosine matrix.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the direction of the axis.
    angle : float
        Angle in radian.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix.
    """
    dcm = get_shiftmat_axis_angle(axis, angle, forward=True)
    return dcm


def any_to_axis_angle(orientation):
    """Converts an orientation to *axis-angle* from any of the following: (1)
    Quaternion, (2) Euler angles, (3) Axis-angle, or (4) Direction cosine
    matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.any_to_quat`.

    Returns
    -------
    axis : (3,) ndarray
        Axis of rotation. This is a unit vector.
    angle : float
        Angle in radian.

    """
    ori_repr = orientation['repr']
    if ori_repr == 'quat':
        quat = np.array(orientation['quat'])
        axis, angle = quat_to_axis_angle(quat)
    elif ori_repr == 'euler':
        euler = np.array(orientation['euler'])
        seq = orientation['seq']
        world = orientation['world']
        axis, angle = euler_to_axis_angle(euler, seq=seq, world=world)
    elif ori_repr == 'axis_angle':
        axis = np.array(orientation['axis'])
        angle = orientation['angle']
    elif ori_repr == 'dcm':
        axis, angle = dcm_to_axis_angle(orientation['dcm'])
    else:
        raise ValueError(
            'Unrecognized orientation repr {0}'.format(ori_repr))
    return axis, angle


def rotate_vector_axis_angle(v, axis, angle):
    """
    Rotates vectors about `axis` by `angle`.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to rotate.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.

    """
    rotmat = get_rotmat_axis_angle(axis, angle)
    return np.dot(v, rotmat.T)


def get_rotmat_axis_angle(axis, angle):
    """Returns the rotation matrix corresponding to an *axis-angle*
    representation.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix.

    """
    R = np.zeros((3,3))
    sin = np.sin(angle)
    cos = np.cos(angle)
    icos = 1.0 - cos
    R[0,0] = axis[0]*axis[0]*icos + cos
    R[0,1] = axis[0]*axis[1]*icos - axis[2]*sin
    R[0,2] = axis[0]*axis[2]*icos + axis[1]*sin
    R[1,0] = axis[0]*axis[1]*icos + axis[2]*sin
    R[1,1] = axis[1]*axis[1]*icos + cos
    R[1,2] = axis[1]*axis[2]*icos - axis[0]*sin
    R[2,0] = axis[2]*axis[0]*icos - axis[1]*sin
    R[2,1] = axis[1]*axis[2]*icos + axis[0]*sin
    R[2,2] = axis[2]*axis[2]*icos + cos
    return R


def extract_axis_angle_from_rotmat(rotmat):
    """
    Extracts axis and angle from a rotation matrix.

    Parameters
    ----------
    rotmat : (3,3) ndarray
        Rotation matrix (must be orthonormal).

    Returns
    -------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian.

    """
    trace = np.trace(rotmat)
    angle = math.acos((trace-1)/2)
    if angle > 0:
        if angle < math.pi:
            u0 = rotmat[2,1] - rotmat[1,2]
            u1 = rotmat[0,2] - rotmat[2,0]
            u2 = rotmat[1,0] - rotmat[0,1]
        else:
            #Find the largest entry in the diagonal of rotmat
            k = np.argmax(np.diag(rotmat))
            if k == 0:
                u0 = math.sqrt(rotmat[0,0]-rotmat[1,1]-rotmat[2,2]+1)/2
                s = 1.0/(2*u0)
                u1 = s*rotmat[0,1]
                u2 = s*rotmat[0,2]
            elif k == 1:
                u1 = math.sqrt(rotmat[1,1]-rotmat[0,0]-rotmat[2,2]+1)/2
                s = 1.0/(2*u1)
                u0 = s*rotmat[0,1]
                u2 = s*rotmat[1,2]
            elif k == 2:
                u2 = math.sqrt(rotmat[2,2]-rotmat[0,0]-rotmat[1,1]+1)/2
                s = 1.0/(2*u2)
                u0 = s*rotmat[0,2]
                u1 = s*rotmat[1,2]
    else:
        u0 = 1.0
        u1 = 0.0
        u2 = 0.0
    return fix_axis_angle(np.array([u0, u1, u2]), angle, normalize=True)


def shift_vector_axis_angle(v, axis, angle, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by an axis-angle representation, shifts vectors from A
    to B or B to A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to shift.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,) or (n,3) ndarray
        Shifted vectors.
    """
    shiftmat = get_shiftmat_axis_angle(axis, angle, forward=forward)
    return np.dot(v, shiftmat.T)


def shift_tensor2_axis_angle(a, axis, angle, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by an axis-angle representation, shifts second order
    tensors from A to B or B to A.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifted second order tensor.
    """
    shiftmat = get_shiftmat_axis_angle(axis, angle, forward=forward)
    return np.einsum('ip,jq,pq', shiftmat, shiftmat, a)


def shift_tensor3_axis_angle(a, axis, angle, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by an axis-angle representation, shifts third order
    tensors from A to B or B to A.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor.
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3,3) ndarray
        Shifted third order tensor.
    """
    shiftmat = get_shiftmat_axis_angle(axis, angle, forward=forward)
    return np.einsum('ip,jq,kr,pqr', shiftmat, shiftmat, shiftmat, a)


def get_shiftmat_axis_angle(axis, angle, forward=False):
    """Returns the shifter matrix corresponding to an *axis-angle*
    representation.

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector along the axis.
    angle : float
        Angle of rotation in radian.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.

    """
    shiftmat = get_rotmat_axis_angle(-axis, angle) 
    if not forward:
        shiftmat = shiftmat.T
    return shiftmat


#DIRECTION COSINE MATRIX-----------------------------------------------
def mat_is_dcm(mat):
    """
    Checks if `mat` is a direction cosine matrix or not.

    Parameters
    ----------
    mat : (3,3) ndarray
        Array to check.

    Returns
    -------
    bool
        ``True`` if `mat` is a direction cosine matrix, ``False`` otherwise.
    """
    return mat_is_rotmat(mat)


def dcm_from_axes(A, B):
    """
    Returns the direction cosine matrix of axes(i.e. frame) B with respect to
    axes(i.e. frame) A.

    Parameters
    ----------
    A : (3,3) ndarray
        The rows of A represent the orthonormal basis vectors of frame A.

    B : (3,3) ndarray
        The rows of B represent the orthonormal basis vectors of frame B.

    Returns
    -------
    (3,3) ndarray
        The dcm of frame B w.r.t. frame A.

    """
    return np.dot(B, A.T)


def dcm_to_quat(dcm):
    """
    Converts a direction cosine matrix to a unit quaternion.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix

    Returns
    -------
    q : (4,) ndarray
        Unit quaternion

    """
    mat = get_rotmat_dcm(dcm)
    axis, angle = extract_axis_angle_from_rotmat(mat)
    return axis_angle_to_quat(axis, angle)


def dcm_to_euler(dcm, seq='XYZ', world=True):
    """
    Converts a direction cosine matrix to a Euler angles.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) ndarray
        The three Euler angles ``phi`` (rotation about X), ``theta`` (rotation
        about Y), and ``psi`` (rotation about Z).
    """
    mat = get_rotmat_dcm(dcm)
    euler = factorize_rotmat(mat, seq=seq, world=world)
    return euler


def dcm_to_axis_angle(dcm):
    """Converts a direction cosine matrix to an *axis-angle* representation.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix

    Returns
    -------
    axis : (3,) ndarray
        Axis of rotation. This is a unit vector.
    angle : float
        Angle in radian.
    """
    mat = get_rotmat_dcm(dcm)
    axis, angle = extract_axis_angle_from_rotmat(mat)
    return (axis, angle)


def any_to_dcm(orientation):
    """Converts an orientation to direction cosine matrix from any of the
    following: (1) Quaternion, (2) Euler angles, (3) Axis-angle, or (4)
    Direction cosine matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.any_to_quat`.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix.
    """
    ori_repr = orientation['repr']
    if ori_repr == 'quat':
        quat = np.array(orientation['quat'])
        dcm = quat_to_dcm(quat)
    elif ori_repr == 'euler':
        euler = np.array(orientation['euler'])
        seq = orientation['seq']
        world = orientation['world']
        dcm = euler_to_dcm(euler, seq=seq, world=world)
    elif ori_repr == 'axis_angle':
        axis = np.array(orientation['axis'])
        angle = orientation['angle']
        dcm = axis_angle_to_dcm(axis, angle)
    elif ori_repr == 'dcm':
        dcm = dcm_to_quat(orientation['dcm'])
    else:
        raise ValueError(
            'Unrecognized orientation repr {0}'.format(ori_repr))
    return dcm


def rotate_vector_dcm(v, dcm):
    """
    Rotates vectors by a direction cosine matrix.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to rotate.
    dcm : (3,3) ndarray
        Direction cosine matrix.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.
    """
    rotmat = get_rotmat_dcm(dcm)
    return np.dot(v, rotmat.T)


def get_rotmat_dcm(dcm):
    """Returns the rotation matrix corresponding to a direction cosine matrix.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix.

    Returns
    -------
    (3,3) ndarray
        Rotation matrix.
    """
    return dcm.T


def shift_vector_dcm(v, dcm, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a direction cosine matrix, shifts vectors from A
    to B or B to A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to shift.
    dcm : (3,3) ndarray
        Direction cosine matrix.   
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,) or (n,3) ndarray
        Shifted vectors.
    """
    shiftmat = get_shiftmat_dcm(dcm, forward=forward)
    return np.dot(v, shiftmat.T)


def shift_tensor2_dcm(a, dcm, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a direction cosine matrix, shifts second-order
    tesnors from A to B or B to A.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor.
    dcm : (3,3) ndarray
        Direction cosine matrix.   
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifted second-order tensor.
    """
    shiftmat = get_shiftmat_dcm(dcm, forward=forward)
    return np.einsum('ip,jq,pq', shiftmat, shiftmat, a)


def shift_tensor3_dcm(a, dcm, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by a direction cosine matrix, shifts third-order
    tesnors from A to B or B to A.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor.
    dcm : (3,3) ndarray
        Direction cosine matrix.   
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3,3) ndarray
        Shifted third-order tensor.
    """
    shiftmat = get_shiftmat_dcm(dcm, forward=forward)
    return np.einsum('ip,jq,kr,pqr', shiftmat, shiftmat, shiftmat, a)


def get_shiftmat_dcm(dcm, forward=False):
    """Returns the shifter matrix corresponding to a direction cosine matrix.

    Parameters
    ----------
    dcm : (3,3) ndarray
        Direction cosine matrix.   
    forward : bool
       Whether to shift forward, i.e., along the orientation or shift reverse.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.
    """
    shiftmat = dcm
    if not forward:
        shiftmat = shiftmat.T
    return shiftmat


#EULER ANGLES-----------------------------------------------------------
def factorize_rotmat(rotmat, seq='XYZ', world=True):
    """Factorize a rotation matrix to obtain the three Euler angles.

    Parameters
    ----------
    rotmat : (3,3) ndarray
        Rotation matrix
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.
    """
    return _factor_rotmat(rotmat, seq=seq, world=world)


def euler_to_euler(euler, seq, world, to_seq, to_world):
    """Convert one set of Euler angles to another.

    Parameters
    ----------
    euler : (3,) ndarray
        The three Euler angles. 
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of `euler`.
    world : bool
        Whether `euler` is with respect to the *world* frame or not.
    to_seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Convert `euler` with sequence `seq` to the sequence `to_seq`.
    to_world : bool
        Whether the converted euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    rotmat = get_rotmat_euler(euler, seq=seq, world=world)
    return factorize_rotmat(rotmat, seq=to_seq, world=to_world)


def euler_to_quat(euler, seq='XYZ', world=True):
    """Convert Euler angles to unit quaternion.

    Parameters
    ----------
    euler : (3,) ndarray
        The three Euler angles. 
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of `euler`.
    world : bool
        Whether `euler` is with respect to the *world* frame or not.

    Returns
    -------
    (4,) ndarray
        Unit quaternion.
    """
    axis, angle = euler_to_axis_angle(euler, seq=seq, world=world)
    return axis_angle_to_quat(axis, angle)


def euler_to_dcm(euler, seq='XYZ', world=True):
    """Convert Euler angles to a direction cosine matrix.

    Parameters
    ----------
    euler : (3,) ndarray
        The three Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of `euler`.
    world : bool
        Whether `euler` is with respect to the *world* frame or not.

    Returns
    -------
    (3,3) ndarray
        Direction cosine matrix
    """
    dcm = get_shiftmat_euler(euler, seq=seq, world=world, forward=True)
    return dcm


def euler_to_axis_angle(euler, seq='XYZ', world=True):
    """Convert Euler angles to an *axis-angle* representation.

    Parameters
    ----------
    euler : (3,) ndarray
        The three Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of `euler`.
    world : bool
        Whether `euler` is with respect to the *world* frame or not.

    Returns
    -------
    axis : (3,) ndarray
        Unit vector along the direction of the axis.
    angle : float
        Angle in radian.    
    """
    rotmat = get_rotmat_euler(euler, seq=seq, world=world)
    axis, angle = extract_axis_angle_from_rotmat(rotmat)
    return (axis, angle)


def any_to_euler(orientation, to_seq, to_world):
    """Converts to Euler angles from any of the following: (1)
    Quaternion, (2) Euler angles, (3) Axis-angle, or (4) Direction cosine
    matrix.

    Parameters
    ----------
    orientation : dict
        See :func:`.any_to_quat`.

    Returns
    -------
    (3,) ndarray
        Euler angles.
    """
    ori_repr = orientation['repr']
    if ori_repr == 'quat':
        quat = np.array(orientation['quat'])
        euler = quat_to_euler(quat, seq=to_seq, world=to_world)
    elif ori_repr == 'euler':
        euler = np.array(orientation['euler'])
        seq = orientation['seq']
        world = orientation['world']
        euler = euler_to_euler(euler, seq, world, to_seq, to_world)
    elif ori_repr == 'axis_angle':
        axis = np.array(orientation['axis'])
        angle = orientation['angle']
        euler = axis_angle_to_euler(axis, angle, seq=to_seq, world=to_world)
    elif ori_repr == 'dcm':
        euler = dcm_to_euler(orientation['dcm'], seq=to_seq, world=to_world)
    else:
        raise ValueError(
            'Unrecognized orientation repr {0}'.format(ori_repr))
    return euler


def rotate_vector_euler(v, euler, seq='XYZ', world=True):
    """
    Rotates vectors with a set of Euler angles.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to rotate.
    euler : (3,) ndarray
        Euler angles.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
        Whether the euler angles are with respect to the *world* frame or not.

    Returns
    -------
    (3,) or (n,3) ndarray
        Rotated vectors.
    """
    rotmat = get_rotmat_euler(euler, seq=seq, world=world)
    return np.dot(v, rotmat.T)


def get_rotmat_euler(euler, seq='XYZ', world=True):
    """
    Returns the rotation matrix for a set of Euler angles.

    Parameters
    ----------
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    """
    return _rotmat_euler(euler, seq=seq, world=world)


def shift_vector_euler(v, euler, seq='XYZ', world=True, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by three Euler angles, shifts vectors from A
    to B or B to A.

    Parameters
    ----------
    v : (3,) or (n,3) ndarray
        A single 3-vector or *n* 3-vectors (the rows of `v`) to shift.
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,) or (n,3) ndarray
        Shifted vectors.
    """
    shiftmat = get_shiftmat_euler(euler, seq=seq, world=world, forward=forward)
    return np.dot(v, shiftmat.T)


def shift_tensor2_euler(a, euler, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by three Euler angles, shifts second-order tensors from A
    to B or B to A.

    Parameters
    ----------
    a : (3,3) ndarray
        A second-order tensor
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifted second order tensor.
    """
    shiftmat = get_shiftmat_euler(euler, forward=forward)
    return np.einsum('ip,jq,pq', shiftmat, shiftmat, a)


def shift_tensor3_euler(a, euler, forward=False):
    """
    Given two frames A and B such that the orientation of frame B with respect
    to frame A is given by three Euler angles, shifts third-order tensors from A
    to B or B to A.

    Parameters
    ----------
    a : (3,3,3) ndarray
        A third-order tensor
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3,3) ndarray
        Shifted third-order tensor.
    """
    shiftmat = get_shiftmat_euler(euler, forward=forward)
    return np.einsum('ip,jq,kr,pqr', shiftmat, shiftmat, shiftmat, a)


def get_shiftmat_euler(euler, seq='XYZ', world=True, forward=False):
    """Returns the shifter matrix for a set of Euler angles.

    Parameters
    ----------
    euler : (3,)
        Euler angles in radian.
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Sequence of Euler angles.
    world : bool
         Whether the euler angles are with respect to the *world* frame or not.
    forward : bool
        If ``True``, shift from A to B. If ``False``, shift from B to A.

    Returns
    -------
    (3,3) ndarray
        Shifter matrix.
    """

    rotmat = get_rotmat_euler(euler, seq=seq, world=world)
    if forward:
        shiftmat = rotmat.T
    else:
        shiftmat = rotmat
    return shiftmat


def _rotmat_euler(euler, seq='XYZ', world=True):
    """
    Returns the rotation matrix for an Euler angle sequence.

    Parameters
    ----------
    euler : (3,) ndarray
        Euler angles
    seq : {'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
        Euler angle sequence.
    world : bool
        Whether the Euler angles are with respect to the *world* frame or not.
    """
    rotmat_funcs = {'XYZ': _rotmat_XYZ, 'XZY': _rotmat_XZY,
                    'YXZ': _rotmat_YXZ, 'YZX': _rotmat_YZX,
                    'ZXY': _rotmat_ZXY, 'ZYX': _rotmat_ZYX
                    }
    if not world:
        euler = -euler
    phi, theta, psi = tuple(euler)
    rotmat = rotmat_funcs[seq](phi, theta, psi)
    if not world:
        rotmat = rotmat.T
    return rotmat


def _factor_rotmat(rotmat, seq='XYZ', world=True):
    factor_rotmat_funcs = {
            'XYZ': _factor_rotmat_XYZ, 'XZY': _factor_rotmat_XZY,
            'YXZ': _factor_rotmat_YXZ, 'YZX': _factor_rotmat_YZX,
            'ZXY': _factor_rotmat_ZXY, 'ZYX': _factor_rotmat_ZYX
            }
    if not world:
        rotmat = rotmat.T
    factors = factor_rotmat_funcs[seq](rotmat) 
    if not world:
        factors = -factors
    return factors


def _rotmat_XYZ(phi, theta, psi):
    rotmat = np.zeros((3,3)) 
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0,0] = cos_theta*cos_psi
    rotmat[0,1] = sin_phi*sin_theta*cos_psi - cos_phi*sin_psi
    rotmat[0,2] = cos_phi*sin_theta*cos_psi + sin_phi*sin_psi
    rotmat[1,0] = cos_theta*sin_psi
    rotmat[1,1] = sin_psi*sin_theta*sin_phi + cos_phi*cos_psi
    rotmat[1,2] = cos_phi*sin_theta*sin_psi - sin_phi*cos_psi
    rotmat[2,0] = -sin_theta
    rotmat[2,1] = sin_phi*cos_theta
    rotmat[2,2] = cos_phi*cos_theta
    return rotmat


def _factor_rotmat_XYZ(rotmat):
    if rotmat[2,0] < 1.0:
        if rotmat[2,0] > -1.0:
            theta = math.asin(-rotmat[2,0])
            psi = math.atan2(rotmat[1,0], rotmat[0,0])
            phi = math.atan2(rotmat[2,1], rotmat[2,2])
        else:
            #Not unique: phi - psi = atan2(-rotmat[1,2], rotmat[1,1])
            theta = math.pi/2
            psi = -math.atan2(-rotmat[1,2], rotmat[1,1])
            phi = 0.0
    else:
        #Not unique: phi + psi = atan2(-rotmat[1,2], rotmat[1,1])
        phi = 0.0
        theta = -math.pi/2
        psi = math.atan2(-rotmat[1,2], rotmat[1,1])
    return np.array([phi, theta, psi])


def _rotmat_XZY(phi, theta, psi):
    rotmat = np.zeros((3,3)) 
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0,0] = cos_theta*cos_psi
    rotmat[0,1] = sin_phi*sin_theta - cos_phi*cos_theta*sin_psi
    rotmat[0,2] = cos_phi*sin_theta + sin_phi*cos_theta*sin_psi
    rotmat[1,0] = sin_psi
    rotmat[1,1] = cos_phi*cos_psi
    rotmat[1,2] = -sin_phi*cos_psi
    rotmat[2,0] = -sin_theta*cos_psi
    rotmat[2,1] = sin_phi*cos_theta + cos_phi*sin_theta*sin_psi
    rotmat[2,2] = cos_phi*cos_theta - sin_phi*sin_theta*sin_psi
    return rotmat


def _factor_rotmat_XZY(rotmat):
    if rotmat[1,0] < 1.0:
        if rotmat[1,0] > -1.0:
            phi = math.atan2(-rotmat[1,2], rotmat[1,1])
            theta = math.atan2(-rotmat[2,0], rotmat[0,0])
            psi = math.asin(rotmat[1,0])
        else:
            #Not unique: phi - theta = atan2(rotmat[2,1], rotmat[2,2])
            phi = 0.0
            theta = -math.atan2(rotmat[2,1], rotmat[2,2])
            psi = -math.pi/2
    else:
        #Not unique: phi + theta = atan2(rotmat[2,1], rotmat[2,2])
        phi = 0.0
        theta = math.atan2(rotmat[2,1], rotmat[2,1])
        psi = math.pi/2
    return np.array([phi, theta, psi])


def _rotmat_YXZ(phi, theta, psi):
    rotmat = np.zeros((3,3)) 
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0,0] = cos_theta*cos_psi - sin_phi*sin_theta*sin_psi
    rotmat[0,1] = -cos_phi*sin_psi
    rotmat[0,2] = sin_theta*cos_psi + sin_phi*cos_theta*sin_psi
    rotmat[1,0] = sin_phi*sin_theta*cos_psi + cos_theta*sin_psi
    rotmat[1,1] = cos_phi*cos_psi
    rotmat[1,2] = sin_theta*sin_psi - sin_phi*cos_theta*cos_psi
    rotmat[2,0] = -cos_phi*sin_theta
    rotmat[2,1] = sin_phi
    rotmat[2,2] = cos_phi*cos_theta
    return rotmat


def _factor_rotmat_YXZ(rotmat):
    if rotmat[2,1] < 1.0:
        if rotmat[2,1] > -1.0:
            phi = math.asin(rotmat[2,1])
            theta = math.atan2(-rotmat[2,0], rotmat[2,2])
            psi = math.atan2(-rotmat[0,1], rotmat[1,1])
        else:
            #Not unique: theta - psi = atan2(rotmat[0,2], rotmat[0,0])
            phi = -math.pi/2
            theta = 0.0
            psi = -math.atan2(rotmat[0,2], rotmat[0,0])
    else:
        #Not unique: theta + psi = atan2(rotmat[0,2], rotmat[0,0])
        phi = math.pi/2
        theta = 0.0
        psi = math.atan2(rotmat[0,2], rotmat[0,0])
    return np.array([phi, theta, psi])


def _rotmat_YZX(phi, theta, psi):
    rotmat = np.zeros((3,3)) 
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0,0] = cos_theta*cos_psi
    rotmat[0,1] = -sin_psi
    rotmat[0,2] = sin_theta*cos_psi
    rotmat[1,0] = sin_phi*sin_theta + cos_phi*cos_theta*sin_psi
    rotmat[1,1] = cos_phi*cos_psi
    rotmat[1,2] = cos_phi*sin_theta*sin_psi - sin_phi*cos_theta
    rotmat[2,0] = sin_phi*cos_theta*sin_psi - cos_phi*sin_theta
    rotmat[2,1] = sin_phi*cos_psi
    rotmat[2,2] = sin_phi*sin_theta*sin_psi + cos_phi*cos_theta
    return rotmat


def _factor_rotmat_YZX(rotmat):
    if rotmat[0,1] < 1.0:
        if rotmat[0,1] > -1.0:
            phi = math.atan2(rotmat[2,1], rotmat[1,1])
            theta = math.atan2(rotmat[0,2], rotmat[0,0])
            psi = math.asin(-rotmat[0,1])
        else:
            #Not unique: theta - phi = atan2(-rotmat[2,0], rotmat[2,2])
            phi = -math.atan2(-rotmat[2,0], rotmat[2,2])
            theta = 0.0
            psi = math.pi/2
    else:
        #Not unique: theta + phi = atan2(-rotmat[2,0], rotmat[2,2])
        phi = math.atan2(-rotmat[2,0], rotmat[2,2])
        theta = 0.0
        psi = -math.pi/2
    return np.array([phi, theta, psi])


def _rotmat_ZXY(phi, theta, psi):
    rotmat = np.zeros((3,3)) 
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0,0] = cos_theta*cos_psi + sin_phi*sin_theta*sin_psi
    rotmat[0,1] = sin_phi*sin_theta*cos_psi - cos_theta*sin_psi
    rotmat[0,2] = cos_phi*sin_theta
    rotmat[1,0] = cos_phi*sin_psi
    rotmat[1,1] = cos_phi*cos_psi
    rotmat[1,2] = -sin_phi
    rotmat[2,0] = sin_phi*cos_theta*sin_psi - sin_theta*cos_psi
    rotmat[2,1] = sin_phi*cos_theta*cos_psi + sin_theta*sin_psi
    rotmat[2,2] = cos_phi*cos_theta
    return rotmat


def _factor_rotmat_ZXY(rotmat):
    if rotmat[1,2] < 1.0:
        if rotmat[1,2] > -1.0:
            phi = math.asin(-rotmat[1,2])
            theta = math.atan2(rotmat[0,2], rotmat[2,2])
            psi = math.atan2(rotmat[1,0], rotmat[1,1])
        else:
            #Not unique: psi - theta = atan2(-rotmat[0,1], rotmat[0,0])
            phi = math.pi/2
            theta = -math.atan2(-rotmat[0,1], rotmat[0,0])
            psi = 0.0
    else:
        #Not unique: psi + theta = atan2(-rotmat[0,1], rotmat[0,0])
        phi = -math.pi/2
        theta = math.atan2(-rotmat[0,1], rotmat[0,0])
        psi = 0.0
    return np.array([phi, theta, psi])


def _rotmat_ZYX(phi, theta, psi):
    rotmat = np.zeros((3,3)) 
    sin_phi = math.sin(phi)
    sin_theta = math.sin(theta)
    sin_psi = math.sin(psi)
    cos_phi = math.cos(phi)
    cos_theta = math.cos(theta)
    cos_psi = math.cos(psi)
    rotmat[0,0] = cos_theta*cos_psi
    rotmat[0,1] = -cos_theta*sin_psi
    rotmat[0,2] = sin_theta
    rotmat[1,0] = sin_phi*sin_theta*cos_psi + cos_phi*sin_psi
    rotmat[1,1] = cos_phi*cos_psi - sin_phi*sin_theta*sin_psi
    rotmat[1,2] = -sin_phi*cos_theta
    rotmat[2,0] = sin_phi*sin_psi - cos_phi*sin_theta*cos_psi
    rotmat[2,1] = sin_phi*cos_psi + cos_phi*sin_theta*sin_psi
    rotmat[2,2] = cos_phi*cos_theta
    return rotmat


def _factor_rotmat_ZYX(rotmat):
    if rotmat[0,2] < 1.0:
        if rotmat[0,2] > -1.0:
            phi = math.atan2(-rotmat[1,2], rotmat[2,2])
            theta = math.asin(rotmat[0,2])
            psi = math.atan2(-rotmat[0,1], rotmat[0,0])
        else:
            #Not unique: psi - phi = atan2(rotmat[1,0], rotmat[1,1])
            phi = -math.atan2(rotmat[1,0], rotmat[1,1])
            theta = -math.pi/2
            psi = 0.0
    else:
        #Not unique: psi + phi = atan2(rotmat[1,0], rotmat[1,1])
        phi = math.atan2(rotmat[1,0], rotmat[1,1])
        theta = math.pi/2
        psi = 0.0
    return np.array([phi, theta, psi])

