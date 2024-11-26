#include <iostream>
#include <math.h>
#include <string>
#include <fstream>

// coded by Jormun 20240922

namespace prs
{
    float pi = 3.1415926, time = 0, dt = 0.001, end = 0.1, ma, re;

    template <int S, int P>
    struct cubit
    {
        cubit *idx[4] = {this, nullptr, nullptr, nullptr};
        int flg[4] = {0, 0, 0, 0};

        float spect[3][4][S]; // r d
        float value[4][P][P], value_x[4][P][P], value_y[4][P][P];
        float base[S][P][P], base_x[S][P][P], base_y[S][P][P];
        float test[S][P][P], test_x[S][P][P], test_y[S][P][P];
        float q[4][2], n[4][2], pos[2][P][P]; // quadrant normal positon

        float legendre(float ps, int k)
        {
            if (k == -1)
                return 0.0;
            else if (k == 0)
                return 1.0;
            else
                return ((2 * k - 1) * ps * legendre(ps, k - 1) - (k - 1) * legendre(ps, k - 2)) / k;
        }

        float degendre(float ps, int k)
        {
            if (k == -1)
                return 0.0;
            else if (k == 0)
                return 0.0;
            else
                return ((2 * k - 1) * (ps * degendre(ps, k - 1) + legendre(ps, k - 1)) - (k - 1) * degendre(ps, k - 2)) / k;
        }

        void clearances(float un[4], int p0, int p1, int flgidx)
        {
            int p0n, p1n;
            if (flg[flgidx] == 0)
                p0n = p0, p1n = p1;
            else if (abs(flg[flgidx]) == 1)
                p0n = 0, p1n = flg[flgidx] > 0 ? p1 : P - 1 - p1;
            else if (abs(flg[flgidx]) == 2)
                p0n = flg[flgidx] > 0 ? p0 : P - 1 - p0, p1n = 0;
            else if (abs(flg[flgidx]) == 3)
                p0n = P - 1, p1n = flg[flgidx] > 0 ? p1 : P - 1 - p1n;

            if (abs(flg[flgidx]) <= 3)
                for (int i = 0; i < 4; ++i)
                    un[i] = idx[flgidx]->value[i][p0n][p1n];

            if (flg[flgidx] == 4)
                ;
        }

        void convection_x(float cx[4], int p0, int p1)
        {
            float u[4], pres, alpha, nn, cxn[4], un[4], presn, alphan;
            u[0] = value[0][p0][p1], u[1] = value[1][p0][p1], u[2] = value[2][p0][p1], u[3] = value[3][p0][p1];
            pres = 0.4 * u[3] - 0.2 * u[1] * u[1] / u[0] - 0.2 * u[2] * u[2] / u[0],
            cx[0] = u[1], cx[1] = u[1] * u[1] / u[0] + pres, cx[2] = u[1] * u[2] / u[0], cx[3] = u[1] / u[0] * (u[3] + pres);
            if (p0 == 0)
                clearances(un, p0, p1, 1), nn = n[1][0];
            else if (p1 == 0)
                clearances(un, p0, p1, 2), nn = n[2][0];
            else if (p0 == P - 1)
                clearances(un, p0, p1, 3), nn = n[3][0];
            else
                return;
            presn = 0.4 * un[3] - 0.2 * un[1] * un[1] / un[0] - 0.2 * un[2] * un[2] / un[0],
            cxn[0] = un[1], cxn[1] = un[1] * un[1] / un[0] + presn, cxn[2] = un[1] * un[2] / un[0], cxn[3] = un[1] / un[0] * (un[3] + presn);
            alpha = sqrt((u[1] * u[1] + u[2] * u[2]) / (u[0] * u[0])) + sqrt(1.4 * pres / u[0]), alphan = sqrt((un[1] * un[1] + un[2] * un[2]) / (un[0] * un[0])) + sqrt(1.4 * presn / un[0]),
            alpha = std::max(alpha, alphan);
            for (int d = 0; d < 4; ++d)
                cx[d] = 0.5 * (cx[d] + cxn[d]) - 0.5 * alpha * nn * (u[d] - un[d]);
        }

        void convection_y(float cy[4], int p0, int p1)
        {
            float u[4], pres, alpha, nn, cyn[4], un[4], presn, alphan;
            u[0] = value[0][p0][p1], u[1] = value[1][p0][p1], u[2] = value[2][p0][p1], u[3] = value[3][p0][p1];
            pres = 0.4 * u[3] - 0.2 * u[1] * u[1] / u[0] - 0.2 * u[2] * u[2] / u[0],
            cy[0] = u[2], cy[1] = u[1] * u[2] / u[0], cy[2] = u[2] * u[2] / u[0] + pres, cy[3] = u[2] / u[0] * (u[3] + pres);
            if (p0 == 0)
                clearances(un, p0, p1, 1), nn = n[1][1];
            else if (p1 == 0)
                clearances(un, p0, p1, 2), nn = n[2][1];
            else if (p0 == P - 1)
                clearances(un, p0, p1, 3), nn = n[3][1];
            else
                return;
            presn = 0.4 * un[3] - 0.2 * un[1] * un[1] / un[0] - 0.2 * un[2] * un[2] / un[0],
            cyn[0] = un[2], cyn[1] = un[1] * un[2] / un[0], cyn[2] = un[2] * un[2] / un[0] + presn, cyn[3] = un[2] / un[0] * (un[3] + presn);
            alpha = sqrt((u[1] * u[1] + u[2] * u[2]) / (u[0] * u[0])) + sqrt(1.4 * pres / u[0]), alphan = sqrt((un[1] * un[1] + un[2] * un[2]) / (un[0] * un[0])) + sqrt(1.4 * presn / un[0]),
            alpha = std::max(alpha, alphan);
            for (int d = 0; d < 4; ++d)
                cy[d] = 0.5 * (cy[d] + cyn[d]) - 0.5 * alpha * nn * (u[d] - un[d]);
        }

        void spect_to_value(int r)
        {
            for (int d = 0; d < 4; ++d)
                for (int p0 = 0; p0 < P; ++p0)
                    for (int p1 = 0; p1 < P; ++p1)
                        value[d][p0][p1] = 0, value_x[d][p0][p1] = 0, value_y[d][p0][p1] = 0;
            for (int d = 0; d < 4; ++d)
                for (int p0 = 0; p0 < P; ++p0)
                    for (int p1 = 0; p1 < P; ++p1)
                        for (int s = 0; s < S; ++s)
                            value[d][p0][p1] += spect[r][d][s] * base[s][p0][p1], value_x[d][p0][p1] += spect[r][d][s] * base_x[s][p0][p1], value_y[d][p0][p1] += spect[r][d][s] * base_y[s][p0][p1];
        }

        void value_to_spect(int r) //
        {
            for (int d = 0; d < 4; ++d)
                for (int s = 0; s < S; ++s)
                    spect[r][d][s] = 0;

            float cx[4], cy[4];

            for (int p0 = 0; p0 < P; ++p0)
                for (int p1 = 0; convection_x(cx, p0, p1), p1 < P - 1; ++p1)
                    for (int d = 0; d < 4; ++d)
                        for (int s = 0; s < S; ++s)
                            spect[r][d][s] += cx[d] * test_x[s][p0][p1];

            for (int p0 = 0; p0 < P; ++p0)
                for (int p1 = 0; convection_y(cy, p0, p1), p1 < P - 1; ++p1)
                    for (int d = 0; d < 4; ++d)
                        for (int s = 0; s < S; ++s)
                            spect[r][d][s] += cy[d] * test_y[s][p0][p1];
        }

        void caculation()
        {
            float l[4], ps[P], ws[P]; // position weight standard
            if (P == 3)
                ps[0] = -1.0, ps[1] = 0.0, ps[2] = +1.0,
                ws[0] = 0.0, ws[1] = 2.0, ws[2] = 0.0;
            else if (P == 4)
                ps[0] = -1.0, ps[1] = -0.57735026918963, ps[2] = +0.57735026918963, ps[3] = +1.0,
                ws[0] = 0.0, ws[1] = 1.0, ws[2] = 1.0, ws[3] = 0.0;
            else if (P == 5)
                ps[0] = -1.0, ps[1] = -0.77459666924148, ps[2] = 0.0, ps[3] = +0.77459666924148, ps[4] = +1.0,
                ws[0] = 0.0, ws[1] = 0.55555555555556, ws[2] = 0.88888888888889, ws[3] = 0.55555555555556, ws[4] = 0.0;
            else if (P == 6)
                ps[0] = -1.0, ps[1] = -0.86113631159405, ps[2] = -0.33998104358486, ps[3] = +0.33998104358486, ps[4] = +0.86113631159405, ps[5] = +1.0,
                ws[0] = 0.0, ws[1] = 0.34785484513745, ws[2] = 0.65214515486255, ws[3] = 0.65214515486255, ws[4] = 0.34785484513745, ws[5] = 0.0;
            else if (P == 7)
                ps[0] = -1.0, ps[1] = -0.90617984593866, ps[2] = -0.53846931010568, ps[3] = 0.0, ps[4] = +0.53846931010568, ps[5] = +0.90617984593866, ps[6] = +1.0,
                ws[0] = 0.0, ws[1] = 0.23692688505619, ws[2] = 0.47862867049937, ws[3] = 0.56888888888889, ws[4] = 0.47862867049937, ws[5] = 0.23692688505619, ws[6] = 0.0;
            else if (P == 18)
                ps[0] = -1.0, ps[1] = -0.9894009349916499, ps[2] = -0.9445750230732326, ps[3] = -0.8656312023878318, ps[4] = -0.7554044083550030, ps[5] = -0.6178762444026438, ps[6] = -0.4580167776572274, ps[7] = -0.2816035507792589, ps[8] = -0.0950125098376374,
                ps[9] = 0.0950125098376374, ps[10] = 0.2816035507792589, ps[11] = 0.4580167776572274, ps[12] = 0.6178762444026438, ps[13] = 0.7554044083550030, ps[14] = 0.8656312023878318, ps[15] = 0.9445750230732326, ps[16] = 0.9894009349916499, ps[17] = 1.0,
                ws[0] = 0.0, ws[1] = 0.0271524594117541, ws[2] = 0.0622535239386479, ws[3] = 0.0951585116824928, ws[4] = 0.1246289712555339, ws[5] = 0.1495959888165767, ws[6] = 0.1691565193950025, ws[7] = 0.1826034150449236, ws[8] = 0.1894506104550685,
                ws[9] = 0.1894506104550685, ws[10] = 0.1826034150449236, ws[11] = 0.1691565193950025, ws[12] = 0.1495959888165767, ws[13] = 0.1246289712555339, ws[14] = 0.0951585116824928, ws[15] = 0.0622535239386479, ws[16] = 0.0271524594117541, ws[17] = 0.0;
            else
                return;
            n[0][0] = q[0][1] - q[1][1], n[0][1] = q[1][0] - q[0][0], n[1][0] = q[1][1] - q[2][1], n[1][1] = q[2][0] - q[1][0],
            n[2][0] = q[2][1] - q[3][1], n[2][1] = q[3][0] - q[2][0], n[3][0] = q[3][1] - q[0][1], n[3][1] = q[0][0] - q[3][0];
            for (int i = 0; i < 4; ++i) //
                l[i] = sqrt(n[i][0] * n[i][0] + n[i][1] * n[i][1]), n[i][0] /= l[i], n[i][1] /= l[i];

            float x0 = (q[0][0] - q[1][0] + q[2][0] - q[3][0]) / 4, x1 = (q[0][0] - q[1][0] - q[2][0] + q[3][0]) / 4, x2 = (q[0][0] + q[1][0] - q[2][0] - q[3][0]) / 4,
                  y0 = (q[0][1] - q[1][1] + q[2][1] - q[3][1]) / 4, y1 = (q[0][1] - q[1][1] - q[2][1] + q[3][1]) / 4, y2 = (q[0][1] + q[1][1] - q[2][1] - q[3][1]) / 4,
                  a = x1 * y0 - x0 * y1, b = x0 * y2 - x2 * y0, c = x1 * y2 - x2 * y1, x_xs, x_ys, y_xs, y_ys, jaco, j2, mass;

            int K = 1, pyramid = 0;
            while ((pyramid += K) < S)
                ++K;
            --K;
            for (int k = 0, s = 0; k <= K; ++k)
                for (int s0 = 0, s1; s1 = k - s0, s0 <= k; ++s, ++s0)
                    for (int p0 = 0; p0 < P; ++p0)
                        for (int p1 = 0; p1 < P; ++p1)
                            x_xs = x0 * ps[p1] + x1, x_ys = x0 * ps[p0] + x2, y_xs = y0 * ps[p1] + y1, y_ys = y0 * ps[p0] + y2,
                            jaco = a * ps[p0] + b * ps[p1] + c, j2 = jaco * jaco, mass = (2.0 * s0 + 1) * (2.0 * s1 + 1) / 4,
                            base[s][p0][p1] = legendre(ps[p0], s0) * legendre(ps[p1], s1) * mass,
                            test[s][p0][p1] = legendre(ps[p0], s0) * legendre(ps[p1], s1) * ws[p0] * ws[p1],
                            base_x[s][p0][p1] = (degendre(ps[p0], s0) * legendre(ps[p1], s1) * y_ys + legendre(ps[p0], s0) * degendre(ps[p1], s1) * (-y_xs)) / jaco * mass,
                            base_y[s][p0][p1] = (degendre(ps[p0], s0) * legendre(ps[p1], s1) * (-x_ys) + legendre(ps[p0], s0) * degendre(ps[p1], s1) * x_xs) / jaco * mass,
                            test_x[s][p0][p1] = (degendre(ps[p0], s0) * legendre(ps[p1], s1) * jaco - legendre(ps[p0], s0) * legendre(ps[p1], s1) * a) / j2 * y_ys * ws[p0] * ws[p1] +
                                                (legendre(ps[p0], s0) * degendre(ps[p1], s1) * jaco - legendre(ps[p0], s0) * legendre(ps[p1], s1) * b) / j2 * (-y_xs) * ws[p0] * ws[p1],
                            test_y[s][p0][p1] = (degendre(ps[p0], s0) * legendre(ps[p1], s1) * jaco - legendre(ps[p0], s0) * legendre(ps[p1], s1) * a) / j2 * (-x_ys) * ws[p0] * ws[p1] +
                                                (legendre(ps[p0], s0) * degendre(ps[p1], s1) * jaco - legendre(ps[p0], s0) * legendre(ps[p1], s1) * b) / j2 * x_xs * ws[p0] * ws[p1];
            for (int k = 0, s = 0; k <= K; ++k)
                for (int s0 = 0, s1; s1 = k - s0, s0 <= k; ++s, ++s0)
                    for (int p = 0; p < P; ++p)
                        jaco = a * ps[p] + b * ps[P - 1] + c,
                        test_x[s][p][P - 1] = legendre(ps[p], s0) * legendre(ps[P - 1], s1) / jaco * l[0] * n[0][0] * ws[p] / 2,
                        test_y[s][p][P - 1] = legendre(ps[p], s0) * legendre(ps[P - 1], s1) / jaco * l[0] * n[0][1] * ws[p] / 2,
                        jaco = a * ps[0] + b * ps[p] + c,
                        test_x[s][0][p] = legendre(ps[0], s0) * legendre(ps[p], s1) / jaco * l[1] * n[1][0] * ws[p] / 2,
                        test_y[s][0][p] = legendre(ps[0], s0) * legendre(ps[p], s1) / jaco * l[1] * n[1][1] * ws[p] / 2,
                        jaco = a * ps[p] + b * ps[0] + c,
                        test_x[s][p][0] = legendre(ps[p], s0) * legendre(ps[0], s1) / jaco * l[2] * n[2][0] * ws[p] / 2,
                        test_y[s][p][0] = legendre(ps[p], s0) * legendre(ps[0], s1) / jaco * l[2] * n[2][1] * ws[p] / 2,
                        jaco = a * ps[P - 1] + b * ps[p] + c,
                        test_x[s][P - 1][p] = legendre(ps[P - 1], s0) * legendre(ps[p], s1) / jaco * l[3] * n[3][0] * ws[p] / 2,
                        test_y[s][P - 1][p] = legendre(ps[P - 1], s0) * legendre(ps[p], s1) / jaco * l[3] * n[3][1] * ws[p] / 2;

            float x, y, xs, ys, xc = 5.0, yc = 5.0, xn, yn, rn, rho, u, v, temp, pres; // eddy value
            for (int p0 = 0; p0 < P; ++p0)
                for (int p1 = 0; p1 < P; ++p1)
                    xs = ps[p0], ys = ps[p1],
                    x = (q[0][0] * (1 + xs) * (1 + ys) + q[1][0] * (1 - xs) * (1 + ys) + q[2][0] * (1 - xs) * (1 - ys) + q[3][0] * (1 + xs) * (1 - ys)) / 4.0,
                    y = (q[0][1] * (1 + xs) * (1 + ys) + q[1][1] * (1 - xs) * (1 + ys) + q[2][1] * (1 - xs) * (1 - ys) + q[3][1] * (1 + xs) * (1 - ys)) / 4.0,
                    xn = x - xc, yn = y - yc, rn = xn * xn + yn * yn, u = 1.0 /*+ 2.5 / pi * exp(0.5 * (1.0 - rn)) * (-yn)*/, v = 1.0 /*+ 2.5 / pi * exp(0.5 * (1.0 - rn)) * (+xn)*/,
                    temp = 1.0 /*- 10.0 / (11.2 * pi * pi) * exp(1.0 - rn)*/, rho = pow(temp, 2.5), pres = pow(rho, 1.4), pos[0][p0][p1] = x, pos[1][p0][p1] = y,
                    value[0][p0][p1] = rho, value[1][p0][p1] = rho * u, value[2][p0][p1] = rho * v, value[3][p0][p1] = 2.5 * pres + 0.5 * rho * (u * u + v * v);

            for (int d = 0; d < 4; ++d)
                for (int s = 0; s < S; ++s)
                    spect[0][d][s] = 0;
            for (int d = 0; d < 4; ++d)
                for (int s = 0; s < S; ++s)
                    for (int p0 = 0; p0 < P; ++p0)
                        for (int p1 = 0; p1 < P; ++p1)
                            spect[0][d][s] += value[d][p0][p1] * test[s][p0][p1];
            spect_to_value(0);
        }
    };

    template <int M, int N, int S, int P>
    struct steps
    {
        cubit<S, P> cbt[N];
        float coord[M][2];
        int compose[N][5];

        steps(std::string file) // ma re
        {
            std::ifstream glasses(file, std::ios::in);
            float trash;
            glasses >> trash;
            for (int m = 0; m < M; ++m)
                glasses >> coord[m][0] >> coord[m][1] >> trash;
            int garbage;
            glasses >> garbage;
            for (int n = 0; n < N; ++n)
                glasses >> compose[n][1] >> compose[n][2] >> compose[n][3] >> garbage, compose[n][0] = compose[n][3], compose[n][4] = compose[n][1]; ///

            for (int n = 0, m; n < N; ++n)
                m = compose[n][1] - 1, cbt[n].q[0][0] = coord[m][0], cbt[n].q[0][1] = coord[m][1], cbt[n].q[1][0] = coord[m][0], cbt[n].q[1][1] = coord[m][1],
                m = compose[n][2] - 1, cbt[n].q[2][0] = coord[m][0], cbt[n].q[2][1] = coord[m][1],
                m = compose[n][3] - 1, cbt[n].q[3][0] = coord[m][0], cbt[n].q[3][1] = coord[m][1], cbt[n].caculation();

            for (int n = 0; n < N; ++n)
                for (int i = 1; i < 4; ++i)
                    if (cbt[n].flg[i] == 0)
                        for (int m = n + 1; m < N; ++m)
                            for (int j = 1; j < 4; ++j)
                                if (cbt[m].flg[j] == 0 && compose[n][i] == compose[m][j + 1] && compose[n][i + 1] == compose[m][j])
                                    cbt[n].idx[i] = &cbt[m], cbt[n].flg[i] = (i == 1 || j == 1) && i != j ? j : -j,
                                    cbt[m].idx[j] = &cbt[n], cbt[m].flg[j] = (i == 1 || j == 1) && i != j ? i : -i;

            for (int n = 0; n < N; ++n)
                for (int i = 1; i < 4; ++i)
                    if (cbt[n].flg[i] == 0)
                        for (int m = n + 1; m < N; ++m)
                            for (int j = 1; j < 4; ++j)
                            {
                                if (cbt[m].flg[j] == 0 && coord[compose[n][i] - 1][0] == coord[compose[n][i + 1] - 1][0])
                                    if (coord[compose[n][i] - 1][1] == coord[compose[m][j + 1] - 1][1] && coord[compose[n][i + 1] - 1][1] == coord[compose[m][j] - 1][1])
                                        cbt[n].idx[i] = &cbt[m], cbt[n].flg[i] = (i == 1 || j == 1) && i != j ? j : -j,
                                        cbt[m].idx[j] = &cbt[n], cbt[m].flg[j] = (i == 1 || j == 1) && i != j ? i : -i;
                                if (cbt[m].flg[j] == 0 && coord[compose[n][i] - 1][1] == coord[compose[n][i + 1] - 1][1])
                                    if (coord[compose[n][i] - 1][0] == coord[compose[m][j + 1] - 1][0] && coord[compose[n][i + 1] - 1][0] == coord[compose[m][j] - 1][0])
                                        cbt[n].idx[i] = &cbt[m], cbt[n].flg[i] = (i == 1 || j == 1) && i != j ? j : -j,
                                        cbt[m].idx[j] = &cbt[n], cbt[m].flg[j] = (i == 1 || j == 1) && i != j ? i : -i;
                            }
        }

        void rk3(float dt)
        {
            for (int n = 0; n < N; ++n)
                cbt[n].value_to_spect(1);
            for (int n = 0; n < N; ++n)
                for (int d = 0; d < 4; ++d)
                    for (int s = 0; s < S; ++s)
                        cbt[n].spect[1][d][s] = cbt[n].spect[0][d][s] + cbt[n].spect[1][d][s] * dt;
            for (int n = 0; n < N; ++n)
                cbt[n].spect_to_value(1);
            for (int n = 0; n < N; ++n)
                cbt[n].value_to_spect(2);
            for (int n = 0; n < N; ++n)
                for (int d = 0; d < 4; ++d)
                    for (int s = 0; s < S; ++s)
                        cbt[n].spect[2][d][s] = (3.0 / 4.0) * cbt[n].spect[0][d][s] + (1.0 / 4.0) * cbt[n].spect[1][d][s] + (1.0 / 4.0) * cbt[n].spect[2][d][s] * dt;
            for (int n = 0; n < N; ++n)
                cbt[n].spect_to_value(2);
            for (int n = 0; n < N; ++n)
                cbt[n].value_to_spect(1);
            for (int n = 0; n < N; ++n)
                for (int d = 0; d < 4; ++d)
                    for (int s = 0; s < S; ++s)
                        cbt[n].spect[0][d][s] = (1.0 / 3.0) * cbt[n].spect[0][d][s] + (2.0 / 3.0) * cbt[n].spect[2][d][s] + (2.0 / 3.0) * cbt[n].spect[1][d][s] * dt;
            for (int n = 0; n < N; ++n)
                cbt[n].spect_to_value(0);
        }

        void vtkview(int interval)
        {
            static int timing = -1;
            ++timing;
            if (timing % interval)
                return;
            std::string file;
            if (timing > 9999)
                file = "view_9999.vtk";
            else if (timing > 999)
                file = "view_" + std::to_string(timing) + ".vtk";
            else if (timing > 99)
                file = "view_0" + std::to_string(timing) + ".vtk";
            else if (timing > 9)
                file = "view_00" + std::to_string(timing) + ".vtk";
            else
                file = "view_000" + std::to_string(timing) + ".vtk";
            std::ofstream pen(file, std::ios::out | std::ios::trunc); //\n?
            pen << "# vtk DataFile Version 2.0" << std::endl
                << "The Penrose Steps are specious." << std::endl
                << "ASCII" << std::endl
                << "DATASET UNSTRUCTURED_GRID" << std::endl
                << "POINTS " << (P - 2) * (P - 2) * N << " float" << std::endl;
            for (int n = 0; n < N; pen << std::endl, ++n)
                for (int p1 = 1; p1 < P - 1; ++p1)
                    for (int p0 = 1; p0 < P - 1; ++p0)
                        pen << cbt[n].pos[0][p0][p1] << ' ' << cbt[n].pos[1][p0][p1] << ' ' << 0 << ' ';
            pen << "CELLS " << (P - 3) * (P - 3) * N << ' ' << (P - 3) * (P - 3) * N * 5 << std::endl;
            for (int n = 0; n < N; ++n)
                for (int p1 = 0; p1 < P - 3; ++p1)
                    for (int p0 = 0; p0 < P - 3; ++p0)
                        pen << 4 << ' ' << n * (P - 2) * (P - 2) + p1 * (P - 2) + p0
                            << ' ' << n * (P - 2) * (P - 2) + p1 * (P - 2) + p0 + 1
                            << ' ' << n * (P - 2) * (P - 2) + (p1 + 1) * (P - 2) + p0 + 1
                            << ' ' << n * (P - 2) * (P - 2) + (p1 + 1) * (P - 2) + p0 << std::endl;
            pen << "CELL_TYPES " << (P - 3) * (P - 3) * N << std::endl;
            for (int n = 0; n < (P - 3) * (P - 3) * N; ++n)
                pen << 9 << ' ';
            pen << std::endl
                << "POINT_DATA " << (P - 2) * (P - 2) * N << std::endl
                << "SCALARS values float 4" << std::endl
                << "LOOKUP_TABLE values" << std::endl;
            for (int n = 0; n < N; pen << std::endl, ++n)
                for (int p1 = 1; p1 < P - 1; ++p1)
                    for (int p0 = 1; p0 < P - 1; ++p0)
                        pen << cbt[n].value[0][p0][p1] << ' ' << cbt[n].value[1][p0][p1] << ' ' << cbt[n].value[2][p0][p1] << ' ' << cbt[n].value[3][p0][p1] << ' ';
        }
    };
}

int main()
{
    prs::steps<100, 162, 6, 18> *stp = new prs::steps<100, 162, 6, 18>("test.dat");
    // prs::steps<33, 48, 6, 18> *stp = new prs::steps<33, 48, 6, 18>("grid.dat");
    stp->vtkview(1);
    while (prs::time < prs::end)
        prs::time += prs::dt, stp->rk3(prs::dt), stp->vtkview(1);
    return 0;
}