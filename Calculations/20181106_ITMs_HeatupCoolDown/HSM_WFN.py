# Module HSM_WFN
"""A module that consists of the methods used to construct the numerical wavefront
aberration from gradients.

"""

from numpy import *
from scipy.sparse import *
from scipy.interpolate import griddata,RectBivariateSpline
from scipy.sparse.linalg import lsmr,lsqr
import gc

def setup_grid(limits,gspacing):
    """Set up square mesh grids for wavefront construction.

    :param limits: the minimum and maximum values of x and y coordinates
    :type limits: *dict*
    :param gspacing: the spacing between two neighboring points on the grid
    :type gspacing: positive *float*
    :Returns: *ndarrays* of x and y grids
    :rtype: *tuple*

    """

    x1,x2 = limits['x']
    y1,y2 = limits['y']

    xr = arange(x1,x2+gspacing,gspacing)
    yr = arange(y1,y2+gspacing,gspacing)

    if xr[-1] > x2:
        xr = xr[:-1]

    if yr[-1] > y2:
        yr = yr[:-1]

    return meshgrid(xr,yr)

def calculate_wf(cents,grads,limits,gspacing,method = 'lsmr'):
    """Calculate the two-dimensional wavefront aberration array from the gradients.

    Note that ``cents`` and ``grads`` must be of the same shape.

    :param cents: the centroids to be used as the positions of the gradients
    :type cents: *ndarray*
    :param grads: the gradients from which the wavefront aberration is to be calcualted
    :type grads: *ndarray*
    :param limits: the minimum and maximum values of x and y coordinates
    :type limits: *dict*
    :param gspacing: the spacing between two neighboring points on the grid
    :type gspacing: positive *float*
    :param method: a method to solve the matrix equation (should be 'lsqr'
                   or 'lsmr')
    :type method: *string*
    :Uses: :func:`~HWS.HSM_WFN.interpolate_gradients`, :func:`~HWS.HSM_WFN.integrate_gradients`
    :Returns: the calculated wavefront
    :rtype: *ndarray*

    """
    gradf_x, gradf_y = interpolate_gradients(cents,grads,limits,gspacing)
    # wf_t = -1.*integrate_grad(gradf_x,gradf_y,gspacing,gspacing)
    wf_t = -1.*integrate_gradients(gradf_x,gradf_y,gspacing,gspacing,method = method)
    wf = wf_t - wf_t.min()

    return wf

def calculate_2nd_grad(cents,grads,limits,gspacing,method = 'lsmr'):
	gradf_x, gradf_y = interpolate_gradients(cents,grads,limits,gspacing)
	print('The shape of grad_x is ' + str(grads.shape))
	print('The shape of gradf_x is ' + str(gradf_x.shape))
	[grad2f_xx,grad2f_xy] = gradient(gradf_x)
	[grad2f_yx,grad2f_yy] = gradient(gradf_y)	
	#vx = zeros_like(gradf_x)
	#vy = zeros_like(gradf_y)
	#for i in range(len(gradf_x)):
	#	for j in range(len(gradf_x)):
	#		vx[i][j] = gradf_x[i][j]/gspacing
	#		vy[i][j] = gradf_y[i][j]/gspacing
	#
	
	return [ [grad2f_xx,grad2f_xy],[grad2f_yx,grad2f_yy]]

def interpolate_gradients(cents,grads,limits,gspacing,method='linear'):
    """Interpolate the gradients.

    This method calls ``griddata`` method in a module ``scipy.interpolate``.

    The input parameter ``method`` must be ``nearest``, ``linear`` or ``cubic``.

    :param cents: the centroids to be used as the positions of the gradients
    :type cents: *ndarray*
    :param grads: the gradients
    :type grads: *ndarray*
    :param limits: the minimum and maximum values of x and y coordinates
    :type limits: *dict*
    :param gspacing: the spacing between two neighboring points on the grid
    :type gspacing: positive *float*
    :param method: a method used to interpolate the gradients
    :type method: *string*
    :Uses: :func:`~HWS.HSM_WFN.setup_grid`
    :Returns: *ndarrays* of x and y components of the interpolated gradients
    :rtype: *tuple*

    """
    mymethod = method
    xg,yg = setup_grid(limits,gspacing)

    gradf_x = griddata(cents,grads[:,0],(xg,yg),method = mymethod)
    gradf_y = griddata(cents,grads[:,1],(xg,yg),method = mymethod)

    gradf_x[isnan(gradf_x)] = 0
    gradf_y[isnan(gradf_y)] = 0

    return (gradf_x,gradf_y)
    
def calculate_wf_prism(limits,gspacing,prism,alpha):
    """Calculate two-dimensional prism wavefront aberration array.

    This method generates a prism wavefront aberration of the magnitude ``prism``
    and the angle ``alpha``.

    :param limits: the minimum and maximum values of x and y coordinates
    :type limits: *dict*
    :param gspacing: the spacing between two neighboring points on the grid
    :type gspacing: positive *float*
    :param prism: the magnitude of the prism
    :type prism: positive *float*
    :param alpha: the angle of the prism
    :type alpha: *float*
    :Uses: :func:`~HWS.HSM_WFN.setup_grid` 
    :Returns: the constructed prism wavefront array
    :rtype: *ndarray*

    """
    xg,yg = setup_grid(limits,gspacing)
    wf = prism * (xg*cos(alpha) + yg*sin(alpha))

def interpolate_wf_default(xc,yc,wf):
    """Interpolate a wavefront aberration array with some default parameters.

    This method can be used to produce a two-dimensional wavefront array from an
    array constructed by :func:`~HWS.HSM_WFN.integrate_gradients`.

    The resulting interpolated wavefront is of a higher resolution (1024 by 1024)
    than the integrated wavefront and thus may be more suitable for a visual
    inspection.

    :param xc: the x coordinates of the wavefront to be interpolated
    :type xc:
    :param yc: the y coordinates of the wavefront to be interpolated
    :type yc:
    :param wf: the wavefront to be interpolated
    :type wf:
    :Uses: :func:`~HS_WFP.setup_grid`
    :Returns: the interpolated wavefron array
    :rtype: *ndarray*
    
    """
    flimits = {}
    flimits['x'] = array([-511.5,511.5])*12e-6
    flimits['y'] = array([-511.5,511.5])*12e-6

    fspacing = 12e-6

    xgf,ygf = setup_grid(flimits,fspacing)

    rbs = RectBivariateSpline(xc,yc,wf.T)
    wf_intp = rbs.ev(xgf,ygf)
    # wf_intp = griddata(coords,wf.flatten(),(xgf,ygf),method = mymethod)

    # if mymethod != 'nearest':
    #     wf_nn = griddata(coords,wf.flatten(),(xg,yg),method = 'nearest')
    #     wf_bi = isnan(wf_intp)
    #     wf_intp[wf_bi] = wf_nn[wf_bi]

    return wf_intp
    
def integrate_gradients(gxi,gyi,dx=1,dy=1,detail=False,method='lsmr'):
    """Numerically integrate the gradients to produce a two-dimensional wavefront
    aberration array.

    :param gxi: the x components of the gradients to be integrated
    :type gxi: *ndarray*
    :param gyi: the y components of the gradients to be integrated
    :type gyi: *ndarray*
    :param dx: the spacing between two nearest points in x direction
    :type dx: positive *float*
    :param dy: the spacing between two nearest points in y direction
    :type dy: positive *float*
    :param detail: a flag to indicate whether to output the integrated wavefront only or
                   to output more information as well
    :param method: the method used to solve the matrix equation (should be either 'lsqr'
                   or 'lsmr')
    :Returns: the integrated wavefront array, (optional) the information on the calculation
    :rtype: *ndarray* or *tuple* of *ndarray* and *dict*

    """
    #  original Matlab code is from:
    #
    #  http://www.mathworks.com/matlabcentral/fileexchange/9734
    #  by John D'Errico
    #  23 Jan 2006 (Updated 27 Jan 2006)
    #  Code covered by the BSD License
    #
    #  Converted to Python by Won Kim
    gx = gxi.copy()
    gy = gyi.copy()
    (ny,nx) = shape(gx)
    # while developing
    # dx = 1
    # dy = 1
    
    dx = tile(dx,[nx-1,1])
    dy = tile(dy,[ny-1,1])
    rhs = zeros(2*nx*ny)
    Af = zeros([2*nx*ny,6])
    L = 0

    indx = 0
    indy = arange(ny)
    ind = indy + indx*ny

    rind = tile(L+arange(ny),[2,1]).T
    cind = array([ind,ind+ny]).T
    dfdx = tile(array([-1., 1.])/dx[0],[ny,1])
    Af[L+arange(ny),:] = concatenate((rind,cind,dfdx),axis = 1)
    rhs[L+arange(ny)] = gx[:,0]
    L = L+ny

    if nx > 2:
        indx,indy = meshgrid(arange(nx-2)+1,arange(ny))
        indx = indx.T.flatten()
        indy = indy.T.flatten()
        ind = indy + indx*ny
        m = ny*(nx-2)

        rind = tile(L+arange(m),[2,1]).T
        cind = array([ind-ny,ind+ny]).T

        dfdx = 1./(dx[indx-1] + dx[indx])
        dfdx = dfdx*array([-1, 1])

        Af[L+arange(m),:] = concatenate((rind,cind,dfdx),axis = 1)
        rhs[L+arange(m)] = gx.T.flatten()[ind]

        L = L + m

    indx = nx - 1
    indy = arange(ny)
    ind = indy + indx*ny
    rind = tile(L+arange(ny),[2,1]).T
    cind = array([ind-ny,ind]).T
    dfdx = tile(array([-1., 1.])/dx[-1],[ny,1])
    Af[L+arange(ny),:] = concatenate((rind,cind,dfdx),axis = 1)
    rhs[L+arange(ny)] = gx[:,-1]
    L = L+ny

    indx = arange(nx)
    indy = 0
    ind = indy + indx*ny
    rind = tile(L+arange(nx),[2,1]).T
    cind = array([ind,ind+1]).T
    dfdy = tile(array([-1,1])/dy[0],[nx,1])
    Af[L+arange(nx),:] = concatenate((rind,cind,dfdy),axis = 1)
    rhs[L+arange(nx)] = gy[1,:]
    L = L + nx

    if ny > 2:
        indx,indy = meshgrid(arange(nx),arange(ny-2)+1)
        indx = indx.flatten()
        indy = indy.flatten()
        ind = indy + indx*ny
        m = nx*(ny-2)

        rind = tile(L+arange(m),[2,1]).T
        cind = array([ind-1,ind+1]).T

        dfdy = 1./(dy[indy-1] + dy[indy])
        dfdy = dfdy*array([-1, 1])

        Af[L+arange(m),:] = concatenate((rind,cind,dfdy),axis = 1)
        rhs[L+arange(m)] = gy.T.flatten()[ind]

        L = L + m
        
    indx = arange(nx)
    indy = ny-1
    ind = indy + indx*ny
    rind = tile(L+arange(nx),[2,1]).T
    cind = array([ind-1,ind]).T
    dfdy = tile(array([-1,1])/dy[-1],[nx,1])
    Af[L+arange(nx),:] = concatenate((rind,cind,dfdy),axis = 1)
    rhs[L+arange(nx)] = gy[-1,:]

    ii = Af[:,0:2].flatten()
    jj = Af[:,2:4].flatten()
    vv = Af[:,4:6].flatten()
    
    A = csc_matrix((vv,(ii,jj)),shape = (2*nx*ny,nx*ny))
    # A = csc_matrix((Af[:,4:6].T,(Af[:,0:2].T,Af[:,2:4].T)),shape = (2*nx*ny,nx*ny))

    # gc.collect()
    # this is where rhs shape changes
    # A0 = A[:,0].copy()
    # dcterm = A0*dc
    #gc.collect()
    # rhs = rhs - dcterm
    # gc.collect()
    # rhs = rhs - dcterm.toarray()
    # rhs = rhs - A.toarray()[:,0]*dc
    # fhat = lsmr(A[:,1:],rhs,atol=0,btol=0,conlim=0)
    gc.collect()
    if method == 'lsqr':
        result = lsqr(A[:,1:],rhs)
    else:
        result = lsmr(A[:,1:],rhs)

    fhat = result[0]
    wf = concatenate((array([0]),fhat)).reshape(ny,nx,order = 'F')

    if detail == True:
        resdat = {}
        if method == 'lsmr':
            resdat['method'] = 'lsmr'
            resdat['x'] = result[0]
            resdat['A'] = A
            resdat['b'] = rhs
            resdat['stop_reason'] = result[1]
            resdat['no_iteration'] = result[2]
            resdat['normr'] = result[3]
            resdat['normar'] = result[4]
            resdat['norama'] = result[5]
            resdat['conda'] = result[6]
            resdat['normx'] = result[7]
        else:
            resdat['method'] = 'lsqr'
            resdat['x'] = result[0]
            resdat['A'] = A
            resdat['b'] = rhs
            resdat['stop_reason'] = result[1]
            resdat['no_iteration'] = result[2]
            resdat['r1norm'] = result[3]
            resdat['r2norm'] = result[4]
            resdat['anorm'] = result[5]
            resdat['acond'] = result[6]
            resdat['arnorm'] = result[7]
            resdat['xnorm'] = result[8]
        return wf,resdat
    else:            
        return wf
