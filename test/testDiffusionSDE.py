# coding=utf-8
import sys

try:
    import unittest
except:
    print("***********************************************************")
    print("* WARNING!! Couldn't import python unittesting framework! *")
    print("* No python tests have been executed                      *")
    print("***********************************************************")
    sys.exit(0)

try:
    import numpy as np
except:
    print("***********************************************************")
    print("* WARNING!! Couldn't import numpy framework! *")
    print("* No python tests have been executed                      *")
    print("***********************************************************")
    sys.exit(-1)

try:
    import crpropa
    from crpropa import nG, kpc, pc, GeV, TeV, PeV, c_light
except Exception as e:
    print("*** CRPropa import failed")
    print(type(e), str(e))
    sys.exit(-2)

def Anderson(x):
    """Scipy independent version of the Anderson-Darling test

    input: x
        - Array of variables to be tested

    output: A2, _Avals, sig
        - A2, test statistic
        - _Avals, list of critical values
        - sig, significance levels in percent

    If the test statistic is larger than the critical value the Null
    hypothesis (The distribution is drawn from normal distribution)
    has to be rejected at the corresponding significane level.

    The implemented version of the Anderson-Darling test is a
    only a test for normality of a given distribution. When it comes
    to testing normality the Anderson-Darling test is one of the most
    powerful tests available. It may be compared to the Kolmogorov-
    Smirnov test but gives more weights to the tails of the distribution.

    Details on the Anderson-Darling test can be found e.g. [1, 2].

    The test is strongly inspired by the scipy [3] implementation
    (see [4]).

    NOTE: The results for the test statistic may differ from the ones
        derived by the scipy implementation. This is due to a different
        treatment of the correction factor. The ratio of critical value
        and test statistic is nevertheless the same.

    [1] https://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test

    [2] http://www.itl.nist.gov/div898/handbook/prc/section2/prc213.htm

    [3] St??fan van der Walt, S. Chris Colbert and Ga??l Varoquaux.
    The NumPy Array: A Structure for Efficient Numerical Computation,
    Computing in Science & Engineering, 13, 22-30 (2011),
    DOI:10.1109/MCSE.2011.37

    [4] https://github.com/scipy/scipy/blob/master/scipy/stats/morestats.py

    [5] Stephens, M A, "EDF Statistics for Goodness of Fit and
    Some Comparisons", Journal of he American Statistical
    Association, Vol. 69, Issue 347, Sept. 1974, pp 730-737

    [6] Merten, L
    The Anderson test statistic was calculated for 6.6e6 independent samples
    of the numpy.random.norm(0,1,N=10000) distribution. The 0.99999-quantile
    of the test statistic distribution is used as the critical value (2.274)
    for the .00001-significance level.
    """


    def CDF_norm(x, mu=0., sigma=1.):

        """Cumulative distribution function

        input: x
            - np.array, list or float

        output: cdf
            - np.array or float

        This function evaluates the cdf of the normal distribution
        at a given position x (see e.g. [1])

        [1] https://en.wikipedia.org/wiki/Normal_distribution
        """

        try:
            x_tilde = (x-mu)/np.sqrt(2.)/sigma
        except TypeError:
            x = np.array(x)
            x_tilde = (x-mu)/np.sqrt(2.)/sigma
        try:
            cdf = [1/2.*(1.+np.math.erf(x)) for x in x_tilde]
            return np.array(cdf)
        except TypeError:
            cdf = 1/2.*(1.+np.math.erf(x_tilde))
            return cdf


    # Critical values from [5, 6]
    _Avals = np.array([0.576, 0.656, 0.787, 0.918, 1.092, 2.274])
    # Significance level in percent
    sig = np.array([15, 10, 5, 2.5, 1, .00001])

    ######################################################

    # The values have to be ordered x1 <= x2 <= ... <= xN
    y = np.sort(x)

    # Rescale the values to a standard normal distribution
    xbar = np.mean(x, axis=0)
    N = len(y)
    s = np.std(x, ddof=1, axis=0)
    w = (y - xbar) / s

    # Calculate the natural logarithm of the cdf and (1-cdf)
    logcdf = np.log(CDF_norm(w))
    logsf = np.log(np.ones(N)-CDF_norm(w))

    # Derive the test statistic
    i = np.arange(1, N + 1)
    A2 = -N - np.sum((2*i - 1.0) / N * (logcdf + logsf[::-1]), axis=0)

    # Rescale it, because mu and sigma were not known
    A2 *= (1 + 4./N - 25./N**2.)



    return {'testStatistic': A2, 'criticalValues': _Avals, 'significanceLevels': sig}

class DiffusionOneDirection(unittest.TestCase):

    ConstMagVec = crpropa.Vector3d(0*nG,0*nG,1*nG)
    BField = crpropa.UniformMagneticField(ConstMagVec)

    precision = 1e-4
    minStep = 1*pc
    maxStep = 10*kpc
    epsilon = 0.

    Dif = crpropa.DiffusionSDE(BField, precision, minStep, maxStep, epsilon)
    DifAdv = crpropa.DiffusionSDE(BField, precision, minStep, maxStep, epsilon)

    def test_Simple(self):
        self.assertEqual(self.Dif.getTolerance(), self.precision)
        self.assertEqual(self.Dif.getMinimumStep(), self.minStep)
        self.assertEqual(self.Dif.getMaximumStep(), self.maxStep)
        self.assertEqual(self.Dif.getEpsilon(), self.epsilon)
        self.assertEqual(self.Dif.getAlpha(), 1./3.) # default Kolmogorov diffusion
        self.assertEqual(self.Dif.getScale(), 1.) # default D(4GeV) = 6.1e28 cm^2/s

    def test_NeutralPropagation(self):
        c = crpropa.Candidate()
        c.current.setId(crpropa.nucleusId(1,0))
        c.current.setEnergy(10*TeV)
        c.current.setDirection(crpropa.Vector3d(1,0,0))

        # check for position
        self.Dif.process(c)
        pos = c.current.getPosition()
        self.assertAlmostEqual(pos.x/pc, 1.) #AlmostEqual due to rounding error
        self.assertEqual(pos.y, 0.)
        self.assertEqual(pos.z, 0.)

        # Step size is increased to maxStep
        self.assertAlmostEqual(c.getNextStep()/self.maxStep, 1.) #AlmostEqual due to rounding error

    def test_NoBFieldPropagation(self):

        ConstMagVec = crpropa.Vector3d(0.)
        BField = crpropa.UniformMagneticField(ConstMagVec)

        precision = 1e-4
        minStep = 1*pc
        maxStep = 10*kpc
        epsilon = 0.

        Dif = crpropa.DiffusionSDE(BField, precision, minStep, maxStep, epsilon)
        c = crpropa.Candidate()
        c.current.setId(crpropa.nucleusId(1,1))
        c.current.setEnergy(10*TeV)
        c.current.setDirection(crpropa.Vector3d(1,0,0))

        # check for position
        Dif.process(c)
        pos = c.current.getPosition()
        self.assertAlmostEqual(pos.x/pc, 1.) #AlmostEqual due to rounding error
        self.assertEqual(pos.y, 0.)
        self.assertEqual(pos.z, 0.)

        # Step size is increased to maxStep
        self.assertAlmostEqual(c.getNextStep()/minStep, 5.) #AlmostEqual due to rounding error


    def test_AdvectivePropagation(self):

        ConstMagVec = crpropa.Vector3d(0.)
        BField = crpropa.UniformMagneticField(ConstMagVec)

        ConstAdvVec = crpropa.Vector3d(0., 1e6, 0.)
        AdvField = crpropa.UniformAdvectionField(ConstAdvVec)

        precision = 1e-4
        minStep = 1*pc
        maxStep = 10*kpc
        epsilon = 0.

        Dif = crpropa.DiffusionSDE(BField, AdvField, precision, minStep, maxStep, epsilon)
        c = crpropa.Candidate()
        c.current.setId(crpropa.nucleusId(1,1))
        c.current.setEnergy(10*TeV)
        c.current.setDirection(crpropa.Vector3d(1,0,0))

        # check for position
        Dif.process(c)
        pos = c.current.getPosition()
        self.assertEqual(pos.x, 0.)
        self.assertAlmostEqual(pos.y, minStep/c_light*1e6)
        self.assertEqual(pos.z, 0.)

        # Step size is increased to maxStep
        self.assertAlmostEqual(c.getNextStep()/minStep, 5.) #AlmostEqual due to rounding error

    msg1 = "Note that this is a statistical test. It might fail by chance with a probabilty O(0.00001)! You should rerun the test to make sure there is a bug."

    def test_DiffusionEnergy10TeV(self):
        x, y, z = [], [], []
        E = 10 * TeV
        D = self.Dif.getScale()*6.1e24*(E/(4*GeV))**self.Dif.getAlpha()
        L_max = 50 * kpc
        std_exp = np.sqrt(2*D*L_max/c_light)
        mean_exp = 0.
        N = 10**4

        maxTra = crpropa.MaximumTrajectoryLength(L_max)

        for i in range(N):
            c = crpropa.Candidate()
            c.current.setId(crpropa.nucleusId(1,1))
            c.current.setEnergy(E)
            while c.getTrajectoryLength() < L_max:
                maxTra.process(c)
                self.Dif.process(c)
                x.append(c.current.getPosition().x)
                y.append(c.current.getPosition().y)
            z.append(c.current.getPosition().z)
        A2 = Anderson(z)['testStatistic']
        self.assertLess(A2, 2.274, msg=self.msg1)
        meanX, meanY, meanZ = np.mean(x), np.mean(y), np.mean(z)
        stdZ = np.std(z)

        # no diffusion in perpendicular direction
        self.assertAlmostEqual(meanX, 0.)
        self.assertAlmostEqual(meanY, 0.)

        # diffusion in parallel direction
        # compare the mean and std of the z-positions with expected values
        # z_mean = 0. (no advection)
        # z_ std = (2*D_parallel*t)^0.5
        # Take 4 sigma errors due to limited number of candidates into account
        stdOfMeans = std_exp/np.sqrt(N)
        stdOfStds = std_exp/np.sqrt(N)/np.sqrt(2.)
        self.assertLess(abs((meanZ-mean_exp)/4./stdOfMeans), 1., msg=self.msg1)
        self.assertLess(abs((stdZ-std_exp)/4./stdOfStds), 1., msg=self.msg1)


    def test_DiffusionEnergy1PeV(self):
        x, y, z = [], [], []
        E = 1 * TeV
        D = self.Dif.getScale()*6.1e24*(E/(4*GeV))**self.Dif.getAlpha()
        L_max = 50 * kpc
        std_exp = np.sqrt(2*D*L_max/c_light)
        mean_exp = 0.
        N = 10**4

        maxTra = crpropa.MaximumTrajectoryLength(L_max)

        for i in range(N):
            c = crpropa.Candidate()
            c.current.setId(crpropa.nucleusId(1,1))
            c.current.setEnergy(E)
            while c.getTrajectoryLength() < L_max:
                maxTra.process(c)
                self.Dif.process(c)
                x.append(c.current.getPosition().x)
                y.append(c.current.getPosition().y)
            z.append(c.current.getPosition().z)
        A2 = Anderson(z)['testStatistic']
        self.assertLess(A2, 2.274, msg=self.msg1)
        meanX, meanY, meanZ = np.mean(x), np.mean(y), np.mean(z)
        stdZ = np.std(z)

        # no diffusion in perpendicular direction
        self.assertAlmostEqual(meanX, 0.)
        self.assertAlmostEqual(meanY, 0.)

        # diffusion in parallel direction
        # compare the mean and std of the z-positions with expected values
        # z_mean = 0. (no advection)
        # z_ std = (2*D_parallel*t)^0.5
        # Take 4 sigma errors due to limited number of candidates into account
        stdOfMeans = std_exp/np.sqrt(N)
        stdOfStds = std_exp/np.sqrt(N)/np.sqrt(2.)
        self.assertLess(abs((meanZ-mean_exp)/4./stdOfMeans), 1., msg=self.msg1)
        self.assertLess(abs((stdZ-std_exp)/4./stdOfStds), 1., msg=self.msg1)

    def test_DiffusionEnergy100PeV(self):
        x, y, z = [], [], []
        E = 10 * PeV
        D = self.Dif.getScale()*6.1e24*(E/(4*GeV))**self.Dif.getAlpha()
        L_max = 50 * kpc
        std_exp = np.sqrt(2*D*L_max/c_light)
        mean_exp = 0.
        N = 10**4

        maxTra = crpropa.MaximumTrajectoryLength(L_max)

        for i in range(N):
            c = crpropa.Candidate()
            c.current.setId(crpropa.nucleusId(1,1))
            c.current.setEnergy(E)
            while c.getTrajectoryLength() < L_max:
                maxTra.process(c)
                self.Dif.process(c)
                x.append(c.current.getPosition().x)
                y.append(c.current.getPosition().y)
            z.append(c.current.getPosition().z)
        A2 = Anderson(z)['testStatistic']
        self.assertLess(A2, 2.274, msg=self.msg1)
        meanX, meanY, meanZ = np.mean(x), np.mean(y), np.mean(z)
        stdZ = np.std(z)

        # no diffusion in perpendicular direction
        self.assertAlmostEqual(meanX, 0.)
        self.assertAlmostEqual(meanY, 0.)

        # diffusion in parallel direction
        # compare the mean and std of the z-positions with expected values
        # z_mean = 0. (no advection)
        # z_ std = (2*D_parallel*t)^0.5
        # Take 4 sigma errors due to limited number of candidates into account
        stdOfMeans = std_exp/np.sqrt(N)
        stdOfStds = std_exp/np.sqrt(N)/np.sqrt(2.)
        self.assertLess(abs((meanZ-mean_exp)/4./stdOfMeans), 1., msg=self.msg1)
        self.assertLess(abs((stdZ-std_exp)/4./stdOfStds), 1., msg=self.msg1)


    def test_FullTransport(self):
        x, y, z = [], [], []
        E = 10 * TeV
        D = self.DifAdv.getScale()*6.1e24*(E/(4*GeV))**self.DifAdv.getAlpha()
        L_max = 50 * kpc
        epsilon = 0.1
        advSpeed = 1e6

        std_exp_x = np.sqrt(2*epsilon*D*L_max/c_light)
        std_exp_y = np.sqrt(2*epsilon*D*L_max/c_light)
        std_exp_z = np.sqrt(2*D*L_max/c_light)
        mean_exp_x = advSpeed * L_max / c_light
        mean_exp_y = 0.
        mean_exp_z = 0.

        N = 10**4

        maxTra = crpropa.MaximumTrajectoryLength(L_max)

        ConstMagVec = crpropa.Vector3d(0*nG,0*nG,1*nG)
        BField = crpropa.UniformMagneticField(ConstMagVec)

        ConstAdvVec = crpropa.Vector3d(advSpeed, 0., 0.)
        AdvField = crpropa.UniformAdvectionField(ConstAdvVec)

        precision = 1e-4
        minStep = 1*pc
        maxStep = 10*kpc

        DifAdv = crpropa.DiffusionSDE(BField, AdvField, precision, minStep, maxStep, epsilon)


        for i in range(N):
            c = crpropa.Candidate()
            c.current.setId(crpropa.nucleusId(1,1))
            c.current.setEnergy(E)
            while c.getTrajectoryLength() < L_max:
                maxTra.process(c)
                DifAdv.process(c)
            x.append(c.current.getPosition().x)
            y.append(c.current.getPosition().y)
            z.append(c.current.getPosition().z)

        # test for normality
        A2 = Anderson(x)['testStatistic']
        self.assertLess(A2, 2.274, msg=self.msg1)
        A2 = Anderson(y)['testStatistic']
        self.assertLess(A2, 2.274, msg=self.msg1)
        A2 = Anderson(z)['testStatistic']
        self.assertLess(A2, 2.274, msg=self.msg1)

        meanX, meanY, meanZ = np.mean(x), np.mean(y), np.mean(z)
        stdX, stdY, stdZ = np.std(x), np.std(y), np.std(z)


        # diffusion in parallel direction
        # compare the mean and std of the z-positions with expected values
        # z_mean = 0. (no advection)
        # z_ std = (2*D_parallel*t)^0.5
        # Take 4 sigma errors due to limited number of candidates into account
        stdOfMeans_x = std_exp_x/np.sqrt(N)
        stdOfStds_x = std_exp_x/np.sqrt(N)/np.sqrt(2.)
        self.assertLess(abs((meanX-mean_exp_x)/4./stdOfMeans_x), 1., msg=self.msg1)
        self.assertLess(abs((stdX-std_exp_x)/4./stdOfStds_x), 1., msg=self.msg1)

        stdOfMeans_y = std_exp_y/np.sqrt(N)
        stdOfStds_y = std_exp_y/np.sqrt(N)/np.sqrt(2.)
        self.assertLess(abs((meanY-mean_exp_y)/4./stdOfMeans_y), 1., msg=self.msg1)
        self.assertLess(abs((stdY-std_exp_y)/4./stdOfStds_y), 1., msg=self.msg1)

        stdOfMeans_z = std_exp_z/np.sqrt(N)
        stdOfStds_z = std_exp_z/np.sqrt(N)/np.sqrt(2.)
        self.assertLess(abs((meanZ-mean_exp_z)/4./stdOfMeans_z), 1., msg=self.msg1)
        self.assertLess(abs((stdZ-std_exp_z)/4./stdOfStds_z), 1., msg=self.msg1)

if __name__ == '__main__':
    unittest.main()

