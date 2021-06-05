from scipy.stats import norm
import copy

def arrange_positions(N, L):
    r = []

    n = round(N**(1/3))
    l = L/(2*n)

    for i in range(n):
        for j in range(n):
            for k in range(n):
                new_atom = [(2*i + 1)*l,  (2*j + 1)*l,  (2*k + 1)*l]
                r.append(new_atom)

    return r

def arrange_velocities(N, T):
    sqrtT = T**0.5
    v = []
    for i in range(N):
        new_atom = [norm.rvs(0, sqrtT) for j in range(3)]
        v.append(new_atom)
    
    return v

def thermostat(T, T_prev, tau):
    return (1 + 1/tau * (T / T_prev  -  1)) ** 0.5

def closest_copy(r, L):
    L_half = L/2

    for i in range(3):
        r_closest_right = r[i] % L

        if r_closest_right <= L_half:
            r[i] = r_closest_right
        else:
            r[i] = r_closest_right - L

def F(r1, r2, L_half, E = None):
    r = [r2[i]-r1[i] for i in range(3)]
    closest_copy(r, L_half)

    r_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]

    r_minus_6 = r_squared**-3
    factor = 24 / r_squared * (r_minus_6 - 2 * r_minus_6*r_minus_6)

    if E is not None:
        E.pot += 4 * (r_minus_6*r_minus_6 - r_minus_6)

    return [r[0] * factor, r[1] * factor,  r[2] * factor]

def accelerations(r, L, E = None):
    N = len(r)
    
    a = [[0]*3 for i in range(N)]

    for i in range(N):
        for j in range(i+1, N):
            a_ij = F(r[i], r[j], L, E)

            a[i][0] += a_ij[0]
            a[i][1] += a_ij[1]
            a[i][2] += a_ij[2]

            a[j][0] -= a_ij[0]
            a[j][1] -= a_ij[1]
            a[j][2] -= a_ij[2]

    return a

def change_positions(r, v, a, dt):
    for i in range(len(r)):
        for j in range(3):
            r[i][j] += v[i][j]*dt  +  0.5 * a[i][j]*dt**2

def temperature(v):
    T = 0
    N = len(v)
    for i in range(N):
        for j in range(3):
            T += v[i][j]**2
    
    T /= (3*N)
    return T


def change_velocities(v, a1, a2, dt, T = -1, tau = 0):
    thermostat_factor = 1

    T_prev = temperature(v)
    if T != -1:
        thermostat_factor = thermostat(T, T_prev, tau)

    for i in range(len(v)):
        for j in range(3):
            v[i][j] += 0.5 * (a1[i][j]  +  a2[i][j]) * dt
            v[i][j] *= thermostat_factor


class Energy:
    pot = 0
    kin = 0

def autocorrelation(v, v_0):
    vv_0_sum = 0
    N = len(v)

    for i in range(N):
        for j in range(3):
            vv_0_sum += v[i][j] * v_0[i][j]
    
    return vv_0_sum / (3 * N)

def delta_r_sqared(r, r_0):
    delta_r2_sum = 0
    N = len(r)

    for i in range(N):
        for j in range(3):
            delta_r2_sum += (r[i][j] - r_0[i][j]) ** 2
    
    return delta_r2_sum / N

####################################

N = 7**3
L = 7.4
T = 1
tau = 10    # thermostat power
tau_0 = 600   # thermostat steps
relaxation = 1000
dt = 0.005 # dt = 0.005 => energy fluctuation ~ 0.5
steps = 2000

xyz_file = open('LJ_animation.xyz', 'w')
exel_file = open('LJ_data.txt', 'w')

r = arrange_positions(N, L)
v = arrange_velocities(N, T)
a = accelerations(r, L)
v_0 = []

for i in range(steps):
    E = Energy()
    if i == relaxation:
        r_0 = copy.deepcopy(r)
    
    if i in (1000, 1200, 1400, 1600, 1800):
        v_0 = copy.deepcopy(v)

    change_positions(r, v, a, dt)
    a_new = accelerations(r, L, E)
    if i < tau_0:
        change_velocities(v, a, a_new, dt, T, tau)
    else:
        change_velocities(v, a, a_new, dt)

    a = a_new

    for j in range(N):
        E.kin += 0.5 * (v[j][0]*v[j][0] + v[j][1]*v[j][1] + v[j][2]*v[j][2])

    # print('{:.10f}, {:.0f}, {:.0f}'.format(E.pot + E.kin,  E.kin, E.pot))

    xyz_file.write(str(N) + '\nframe ' + str(i) + '\n')
    for j in range(N):
        xyz_file.write('a ' + ' '.join(map(str, r[j])) + '\n')
    
    exel_file.write(str(i))
    exel_file.write(' ')
    exel_file.write(str(temperature(v)))
    exel_file.write(' ')
    if i >= relaxation:
        exel_file.write(str(autocorrelation(v, v_0)))
        exel_file.write(' ')
        exel_file.write(str(delta_r_sqared(r, r_0)))
    else:
        exel_file.write('0 0')
    exel_file.write('\n')
    
    if i%100 == 0:
        print(i)