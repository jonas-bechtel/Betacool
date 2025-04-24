#include "xForce.h"

xFrParam::xFrParam()
{
    V_tr_e._(m_ / s_);
    // ----- 05.06  for 3D force
    V_tr_x._(m_ / s_);
    V_tr_y._(m_ / s_);
    //-----
    V_long_e._(m_ / s_);
    Ttemp._(eV_);
    Ttemp_centre._(eV_);
    Ltemp._(eV_);
    TempEff._(eV_);
    Theta_Eff._(U1_);
    V_eff_e._(m_ / s_);
    Smoos._(U1_);
    mfield._(k_3 * G_);
    n_e._(1 / (m_ ^ 3));
    tau._(s_);
    T_plasma._(s_);
    //Undulator parameters
    lambda._(cm_);
    B_field._(G_);
    r_0._(u_6 * m_);
    Theta_U._(U1_);
    V_und._(m_ / s_);
}

xForce::xForce()
{
   type = 0; //type of friction force formula
   delta._((Q_ ^ 2) / (s_ ^ 2));
   Vtr_min._(m_ / s_);
   Vtr_max._(m_ / s_);
   Vlong_min._(m_ / s_);
   Vlong_max._(m_ / s_);
   V_max._(m_ / s_);
   V_min._(m_ / s_);
   Velocity._(m_ / s_);

   Vtr._(m_ / s_);
   Ftr._(eV_ / m_);
   for (int i = 0; i < 3; i++)
   {
      v[i]._(m_ / s_);
      f[i]._(eV_ / m_);
   }
}

void xForce::Budker(xFrParam F)
{
   delta = 4 * U_pi * F.n_e * F.Z * F.Z * (U_e ^ 4) / U_me; //Friction constant
   doubleU R_max(m_);                                       // max impact parameter
   doubleU R_max1(m_);                                      // for max impact parameter (comparisom)
   doubleU R_max2(m_);
   doubleU R_min(m_);                    // min impact parameter
   doubleU E_e(M_6 * eV_ * 0.5110034);   // rest electron energy
   doubleU theta(m_ / s_);               // total ion velocity
   doubleU L_C;                          // Coloumb logarithms
   doubleU delta_x((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_s((s_ ^ 2) / (m_ ^ 2)); //
   double x;

   theta = ((Vtr ^ 2) + (v[2] ^ 2)) ^ 0.5;

   if (theta() > 0.)
   {
      //components of the friction force
      delta_x = 0;
      delta_s = 0;
      //maximum and minimum impact parameters constant
      R_max = (((theta * theta + F.V_tr_e * F.V_tr_e) * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);
      R_max1 = ((3 * F.Z / F.n_e) ^ (1. / 3.));
      if (R_max < R_max1)
         R_max = R_max1;
      R_max2 = theta * F.tau;

      if (R_max() > R_max2())
         R_max = R_max2;

      R_min = F.Z * U_e * U_e / (U_me * (theta * theta + F.V_tr_e * F.V_tr_e));

      //       	R_min /= (theta + F.V_tr_e)*(theta + F.V_tr_e);
      if (R_max() > R_min())
         L_C = U_Ln(R_max / R_min);
      else
         L_C = 0;

      x = theta() / F.V_tr_e();
      delta_x =  Vtr * L_C * phi(x) / (theta * theta * theta);
      delta_s = v[2] * L_C * phi(x) / (theta * theta * theta);

      Ftr  = -delta_x * delta;
      f[2] = -delta_s * delta;
   }
   else 
   { // total angle is 0
      Ftr  = 0;
      f[2] = 0; 
   }
}

void xForce::NonMag(xFrParam F)
{
   delta = 4 * U_pi * F.n_e * F.Z * F.Z * (U_e ^ 4) / U_me; //Friction constant
   doubleU R_max(m_);                                       // max impact parameter
   doubleU R_max1(m_);                                      // for max impact parameter (comparisom)
   doubleU R_max2(m_);
   doubleU R_max3(m_);
   doubleU R_min(m_);                  // min impact parameter
   doubleU E_e(M_6 * eV_ * 0.5110034); // rest electron energy
   doubleU theta(m_ / s_);             // total ion velocity
   doubleU theta_e(m_ / s_);
   doubleU norm(m_ / s_);
   doubleU L_C;  // Coloumb logarithms
   doubleU L_C3; // Coloumb logarithms
   doubleU Maxw;
   //--------------------------------------------29.11.05------------------
   doubleU B_long;
   doubleU B_trans;
   doubleU DelExp;
   doubleU B_long_step;
   doubleU B_trans_step;
   doubleU q;
   doubleU Vlong_sig;
   doubleU Vtr_sig;
   doubleU sig_sig;
   //--------------------------------------------29.11.05------------------
   doubleU delta_x((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_s((s_ ^ 2) / (m_ ^ 2)); //
   doubleU d_x((s_) / (m_));             //
   doubleU d_s((s_) / (m_));             //
   int i, j, l;
   //int nfi = 15;
   double dfi;
   double fi0;
   dfi = M_PI / nfi;
   //dfi = 2.0*M_PI/nfi;
   fi0 = dfi / 2;
   doubleU v_l(m_ / s_);               // step over longitudinal electron velocity
   doubleU v_t(m_ / s_);               // step over transverse electron velocity
   doubleU U_rel((m_ ^ 2) / (s_ ^ 2)); // square of relative velocity

   theta = ((Vtr ^ 2) + (v[2] ^ 2)) ^ 0.5;
   theta_e = ((F.V_tr_e ^ 2) + (F.V_long_e ^ 2)) ^ 0.5;

   if (theta() > 0.)
   {
      //components of the friction force
      delta_x = 0;
      delta_s = 0;

      //maximum and minimum impact parameters constant
      R_max1 = (3 * F.Z / F.n_e) ^ (1. / 3.);
      R_max2 = theta * F.tau;
      R_max3 = (F.V_tr_e * F.V_tr_e * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5;
      if (smooth)
         R_max = ((theta * theta + theta_e * theta_e) * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5;
      else
      {
         if (theta > theta_e)
            R_max = (theta * theta * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5;
         else
            R_max = (theta_e * theta_e * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5;
      }
      if (R_max < R_max1)
         R_max = R_max1; // choose max
      if (R_max > R_max2)
         R_max = R_max2; // choose min

      if (asimptotic) // Asimptotic -------------------------------------------
      {
         R_min = F.Z * U_e * U_e / (U_me * (theta * theta + theta_e * theta_e));

         if (F.undulator && R_min < F.r_0)
            R_min = F.r_0;

         if (R_max > R_min)
            L_C = U_Ln(R_max / R_min);
         else
            L_C = 0;

         if (R_max3 > R_min)
            L_C3 = U_Ln(R_max3 / R_min);
         else
            L_C3 = 0;

         if (theta >= F.V_tr_e) // I //
         {
            delta_x = Vtr * L_C / (theta * theta * theta);
            delta_s = v[2] * L_C / (theta * theta * theta);
            if ((asimptotic == 1) && (v[2]()))
               delta_s -= L_C3 / (v[2] * v[2]) * sqrt(2. / M_PI);
         }
         else if (theta > F.V_long_e) // II //
         {
            delta_x = Vtr * L_C / (F.V_tr_e * F.V_tr_e * F.V_tr_e);
            if (v[2]() == 0.0)
               delta_s = 0.;
            else
            {
               if (asimptotic == 2)
                  delta_s = ((v[2]() < 0 ? -1 : 1) * L_C / (F.V_tr_e * F.V_tr_e));
               else
                  delta_s = (v[2] * L_C / (((v[2] * v[2] + F.V_long_e * F.V_long_e) ^ 0.5) * F.V_tr_e * F.V_tr_e)) - (v[2] * L_C3 / (F.V_tr_e * F.V_tr_e * F.V_tr_e) * sqrt(2. / M_PI));
            }
         }
         else
         {
            if (asimptotic == 2) // III //
               delta_s = (v[2] * L_C) / (F.V_long_e * F.V_tr_e * F.V_tr_e);
            else
               delta_s = (v[2] * L_C) / (((v[2] * v[2] + F.V_long_e * F.V_long_e) ^ 0.5) * F.V_tr_e * F.V_tr_e);
         }
      }
      else // Numerical --------------------------------------------
      {
         //--------------------------------------------29.11.05------------------
         if (DelPopolo)
         {
            delta = 2 * ((2 * U_pi) ^ 0.5) * F.n_e * F.Z * F.Z * (U_e ^ 4) / U_me;
            R_min = F.Z * U_e * U_e / (U_me * (theta * theta + theta_e * theta_e));
            if (F.undulator && R_min < F.r_0)
               R_min = F.r_0;
            if (R_max() > R_min())
               L_C = U_Ln(R_max / R_min);
            else
               L_C = 0;
            Vlong_sig = (v[2] * v[2]) / (2 * F.V_tr_e * F.V_tr_e);
            Vtr_sig = (Vtr * Vtr) / (2 * F.V_tr_e * F.V_tr_e);
            sig_sig = (F.V_long_e * F.V_long_e) / (F.V_tr_e * F.V_tr_e);
            q = 0.0;
            B_trans = 0.0;
            B_long = 0.0;
            do
            {
               // --- 05.06
               if ((Vtr_sig / (1.0 + q)) + (Vlong_sig / (sig_sig + q))(U1_) < 300)
               {
                  DelExp = U_Exp((Vtr_sig / (1.0 + q)) + (Vlong_sig / (sig_sig + q)));
                  //B_long_step = 1.0/(DelExp*(1.0+q)*(1.0+q)*((sig_sig+q)^1.5));
                  B_long_step = 1.0 / (DelExp * (1.0 + q) * ((sig_sig + q) ^ 1.5));
                  B_trans_step = 1.0 / (DelExp * (1.0 + q) * (1.0 + q) * ((sig_sig + q) ^ 0.5));
               }
               else
               {
                  B_long_step = 0;
                  B_trans_step = 0;
               }
               B_trans += B_trans_step * q_step;
               B_long += B_long_step * q_step;
               q += q_step;
            } while (q() < q_max);

            delta_x = Vtr * B_trans * L_C / (F.V_tr_e * F.V_tr_e * F.V_tr_e);
            delta_s = v[2] * B_long * L_C / (F.V_tr_e * F.V_tr_e * F.V_tr_e);
         }
         else
         {
            v_t = (F.V_tr_e * 3) / dt;
            v_l = (F.V_long_e * 3) / dl;
            d_x = 0;
            d_s = 0;
            norm = 0;

            for (i = 1 - dl; i < dl; i++)
               for (j = 0; j < dt; j++)
                  for (l = 0; l < nfi; l++)
                  {
                     U_rel = (((v[2] - (v_l * i)) * (v[2] - (v_l * i)) +
                               (Vtr - ((v_t * j) * cos(dfi * l + fi0))) * (Vtr - ((v_t * j) * cos(dfi * l + fi0))) +
                               (v_t * j) * sin(dfi * l + fi0) * (v_t * j) * sin(dfi * l + fi0)));

                     if (U_rel() != 0.0)
                     {
                        if (rhomin)
                           R_min = F.Z * U_e * U_e / (U_me * U_rel);
                        else
                        {
                           if (rmsplus)
                              R_min = F.Z * U_e * U_e / (U_me * (theta * theta + theta_e * theta_e));
                           else
                              R_min = F.Z * U_e * U_e / (U_me * theta * theta);
                        }
                        if (F.undulator && R_min < F.r_0)
                           R_min = F.r_0;
                        if (R_max() > R_min())
                           L_C = U_Ln(R_max / R_min);
                        else
                           L_C = 0;

                        Maxw = U_Exp(((v_t * j) * (v_t * j)) / (F.V_tr_e * F.V_tr_e * 2.) + ((v_l * i) * (v_l * i)) / (F.V_long_e * F.V_long_e * 2.));
                        d_x += (Vtr - ((v_t * j) * cos(dfi * l + fi0))) * L_C * (v_t * j) / ((U_rel ^ 1.5) * Maxw);
                        d_s += (v[2] - (v_l * i)) * L_C * (v_t * j) / ((U_rel ^ 1.5) * Maxw);
                        norm += (v_t * j) / Maxw;
                     }
                  }
            delta_x = d_x / norm;
            delta_s = d_s / norm;
           
         }
      }
      //--------------------------------------------29.11.05------------------
      Ftr = -delta_x * delta;
      f[2] = -delta_s * delta;
   }
   else
   {
      Ftr = 0;
      f[2] = 0; // total angle is 0
   }
}

void xForce::DerSkr(xFrParam F)
{
   delta = 2 * U_pi * F.n_e * F.Z * F.Z * (U_e ^ 4) / U_me; //Friction constant
   double t;
   doubleU R_max(m_);  // max impact parameter
   doubleU R_max1(m_); // for max impact parameter (comparisom)
   doubleU R_max2(m_);
   doubleU R_min(m_); // min impact parameter
   doubleU R_f(m_);   // intermediate impact parameter
   doubleU ro_e(m_);  // electron Larmour radius
   //     05.06
   doubleU R_min_mag(m_);              // min impact for magnetized collisions
   doubleU E_e(M_6 * eV_ * 0.5110034); // rest electron energy
   doubleU k_xz, k_s;                  // koefficients Kx,z
   doubleU theta_e(m_ / s_);           // total electron velocity spread
   doubleU theta(m_ / s_);             // total ion velocity
   doubleU L_M, L_A, L_F;              // Coloumb logarithms
   doubleU N_col;                      // number of collisions
   // Numerical integration
   doubleU delta_x_m((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_s_m((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_x_a((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_s_a((s_ ^ 2) / (m_ ^ 2)); //
   doubleU Maxw;
   doubleU delta_x((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_s((s_ ^ 2) / (m_ ^ 2)); //
   doubleU d_x((s_) / (m_));             //
   doubleU d_s((s_) / (m_));             //
   doubleU d_x_m((s_ ^ 2) / (m_ ^ 2));   //
   doubleU d_s_m((s_ ^ 2) / (m_ ^ 2));   //
   doubleU norm_m;
   int i, j, l, n;
   int nfi = 7;
   double dfi;
   double fi0;
   dfi = M_PI / nfi;
   fi0 = dfi / 2;
   double y, z;
   double alpha, sum, sum2;
   doubleU v_l(m_ / s_);               // step over longitudinal electron velocity
   doubleU v_e(m_ / s_);               // longitudinal electron velocity
   doubleU v_t(m_ / s_);               // step over transverse electron velocity
   doubleU U_rel((m_ ^ 2) / (s_ ^ 2)); // square of relative velocity
   doubleU norm(m_ / s_);
   doubleU ellipsoid; // constant determined region of velocity

   //Components of the ion velocity
   theta = ((Vtr ^ 2) + (v[2] ^ 2)) ^ 0.5;

   if (theta() > 0.)
   {
      //components of the friction force
      delta_x = 0;
      delta_s = 0;
      delta_x_m = 0;
      delta_s_m = 0;
      delta_x_a = 0;
      delta_s_a = 0;
      //  05.06
      //maximum and minimum impact parameters constant
      R_max = (((theta * theta + F.V_long_e * F.V_long_e) * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);
      R_max1 = ((3 * F.Z / F.n_e) ^ (1. / 3.));
      R_max2 = ((theta * theta + F.V_long_e * F.V_long_e) ^ 0.5) * F.tau;

      //If Debue sphere does not contain enough electrons R_max is determined by electron density
      if (R_max() < R_max1())
         R_max = R_max1;
      if (R_max() > R_max2())
         R_max = R_max2;
      // 05.06
      R_min = F.Z * U_e * U_e / (U_me * (theta * theta + F.V_long_e * F.V_long_e));
      //R_min = F.Z*U_e*U_e/(U_me*(theta * theta + F.V_tr_e*F.V_tr_e));
      k_xz = 1 - 3 * (v[2] / theta) * (v[2] / theta);
      k_s = 2 + k_xz;

      theta_e = ((F.V_tr_e * F.V_tr_e + F.V_long_e * F.V_long_e) ^ 0.5);
      ro_e = U_me * F.V_tr_e * U_c / (U_e * F.mfield);
      R_f = ro_e * ((theta * theta + F.V_long_e * F.V_long_e) ^ 0.5) / F.V_tr_e;
      //          05.06
      R_min_mag = F.Smoos * ro_e;
      if (R_min_mag() < R_min())
         R_min_mag = R_min;
      L_M = 0;
      L_A = 0;
      L_F = 0;
      if (R_max > R_min_mag)
         L_M = U_Ln(R_max / R_min_mag);
      if ((F.Smoos * ro_e) > R_f)
         L_A = U_Ln(F.Smoos * ro_e / R_f);
      if (R_f > R_min)
         L_F = U_Ln(R_f / R_min);

      N_col = 1 + F.V_tr_e / (U_pi * ((theta * theta + F.V_long_e * F.V_long_e) ^ 0.5));

      if (Mag_As) // Asimptotic ----------------------------------------------
      {
         if (theta >= F.V_tr_e) /* I */
         {
            delta_x = Vtr * (2 * L_F + k_xz * L_M) / (theta ^ 3);
            delta_s = v[2] * (2 * L_F + k_s * L_M + 2) / (theta ^ 3);
         }
         else
         {
            if (theta <= F.V_long_e) /* III */
            {
               if (Vtr() == 0)
                  delta_x = 0;
               else
                  delta_x = Vtr * (2 * (L_F + N_col * L_A) / (F.V_tr_e ^ 3) + (L_M * U_Ln(F.V_long_e / Vtr)) / (F.V_long_e ^ 3));
               delta_s = v[2] * (2 * (L_F + N_col * L_A) / (F.V_tr_e ^ 2) + (L_M / (F.V_long_e ^ 2))) / F.V_long_e;
            }
            else /* II_a + II_b */
            {
               delta_x = Vtr * ((2 * (L_F + N_col * L_A)) / (F.V_tr_e ^ 3) + (k_xz * L_M / (theta ^ 3)));
               ellipsoid = (Vtr ^ 2) / (F.V_tr_e ^ 2) + ((v[2] ^ 2) / (F.V_long_e ^ 2));
               if (ellipsoid <= 1.) /* II_b */
                  delta_s = v[2] * ((2 * (L_F + N_col * L_A)) / (F.V_tr_e ^ 2) / F.V_long_e) + ((v[2] * (k_s * L_M + 2)) / (theta ^ 3));
               else /* II_a */
               {    /* Signum */
                  t = v[2]() < 0. ? -1. : 1.;
                  delta_s = (2 * t * (L_F + N_col * L_A)) / (F.V_tr_e ^ 2) + ((v[2] * (k_s * L_M + 2)) / (theta ^ 3));
                  if (v[2]() == 0.0)
                     delta_s = 0.;
               }
            }
         }
      }
      else // Numerical -----------------------------------------------
      {
         if (Fast)
         { // Non-magnitized -------------------------------------------
            if (R_f > R_min)
            {
               L_F = U_Ln(R_f / R_min);
               if (R_f > R_max)
                  L_F = U_Ln(R_max / R_min);

               v_t = (F.V_tr_e * 3) / dt;
               v_l = (F.V_long_e * 3) / dl;
               n = 0;
               d_x = 0;
               d_s = 0;
               norm = 0;

               for (i = 1 - dl; i < dl; i++)
                  for (j = 0; j < dt; j++)
                     for (l = 0; l < nfi; l++)
                     {
                        U_rel = (((v[2] - (v_l * i)) * (v[2] - (v_l * i)) +
                                  (Vtr - ((v_t * j) * cos(dfi * l + fi0))) * (Vtr - ((v_t * j) * cos(dfi * l + fi0))) +
                                  (v_t * j) * sin(dfi * l + fi0) * (v_t * j) * sin(dfi * l + fi0)));

                        if (U_rel() != 0.0)
                        {
                           R_min = F.Z * U_e * U_e / (U_me * U_rel);
                           if (R_min < R_f)
                           {
                              L_F = U_Ln(R_f / R_min);
                              if (R_f > R_max)
                                 L_F = U_Ln(R_max / R_min);
                           }
                           else
                              L_F = 0.0;

                           Maxw = U_Exp(((v_t * j) * (v_t * j)) / (F.V_tr_e * F.V_tr_e * 2.) + ((v_l * i) * (v_l * i)) / (F.V_long_e * F.V_long_e * 2.));
                           d_x += (Vtr - ((v_t * j) * cos(dfi * l + fi0))) * L_F * (v_t * j) / ((U_rel ^ 1.5) * Maxw);
                           d_s += (v[2] - (v_l * i)) * L_F * (v_t * j) / ((U_rel ^ 1.5) * Maxw);
                           norm += (v_t * j) / Maxw;
                           n++;
                        }
                     }
               delta_x = d_x / norm * 2.0;
               delta_s = d_s / norm * 2.0;
            }
         }

         if (Magnetized)
         {
            if (R_min_mag > R_max)
               L_M = 0.0;
            else
               L_M = U_Ln(R_max / R_min_mag);
            //L_M = U_Ln((R_max+R_min_mag)/R_min_mag);
            if (Pestrikov)
            { // Pestrikov -------------------------------------------
               z = fabs(Vtr() / F.V_long_e());
               y = v[2]() / F.V_long_e();
               dfi = M_PI / nfiP;
               fi0 = dfi / 2;
               alpha = fi0 - M_PI / 2.0;
               sum = 0.0;
               sum2 = 0.0;
               for (l = 0; l < nfiP; l++)
               {
                  sum += (y * cos(alpha) + z * sin(alpha)) * exp(-0.5 * (y + z * tan(alpha)) * (y + z * tan(alpha)));
                  sum2 += tan(alpha) * (y * cos(alpha) + z * sin(alpha)) * exp(-0.5 * (y + z * tan(alpha)) * (y + z * tan(alpha)));
                  alpha += dfi;
               }
               sum *= dfi;
               sum2 *= dfi;
               delta_s_m = L_M * sum / (((2.0 * U_pi) ^ 0.5) * F.V_long_e * F.V_long_e);
               delta_x_m = L_M * sum2 / (((2.0 * U_pi) ^ 0.5) * F.V_long_e * F.V_long_e);
               if (Vtr() <= 0.0)
                  delta_x_m *= -1.0;
            }
            else if (theta <= 0.003 * F.V_long_e)
            {
               if (Vtr() != 0.0)
                  delta_x_m = 2.0 * L_M * Vtr * U_Ln(F.V_long_e / U_Abs(Vtr)) / (((2.0 * U_pi) ^ 0.5) * F.V_long_e * F.V_long_e * F.V_long_e);
               delta_s_m = 2.0 * L_M * v[2] / (((2.0 * U_pi) ^ 0.5) * F.V_long_e * F.V_long_e * F.V_long_e);
            }
            else
            { // Magnetized -------------------------------------------
               v_l = (F.V_long_e * 4.1) / N_M;
               n = 0;
               d_x_m = 0;
               d_s_m = 0;
               norm_m = 0;
               v_e = v[2] + (v_l / 2.0);

               do
               { //L_M = U_Ln((R_max+R_min_mag)/R_min_mag);
                  U_rel = (((v[2] - v_e) * (v[2] - v_e) + Vtr * Vtr));
                  Maxw = U_Exp(-1.0 * (v_e * v_e) / (F.V_long_e * F.V_long_e * 2.));
                  norm_m += Maxw;
                  /*
                  R_min = F.Z*U_e*U_e/(U_me*(U_rel+F.V_tr_e*F.V_tr_e));
                  if (R_min > R_max) L_M = 0.0;
                  else
                  {  if (R_min > ro_e)L_M = U_Ln((R_max+R_min)/R_min);
                  }
                  */
                  if (R_min_mag > R_max)
                     L_M = 0.0;
                  else
                  {
                     L_M = U_Ln(R_max / R_min_mag);
                     //if(R_min > ro_e)L_M = U_Ln((R_max+R_min)/R_min);
                  }
                  d_x_m += Vtr * (Vtr * Vtr - 2.0 * (v[2] - v_e) * (v[2] - v_e)) * L_M * Maxw / (U_rel ^ 2.5);
                  d_s_m += 3.0 * Vtr * Vtr * (v[2] - v_e) * L_M * Maxw / (U_rel * (U_rel ^ 1.5));

                  if (Constant)
                  {
                     if ((U_rel ^ 0.5) > (1.0 * F.V_long_e))
                        d_s_m += 2.0 * (v[2] - v_e) * Maxw / (U_rel ^ 1.5);
                  }
                  v_e += v_l;
               } while (v_e <= (4.1 * F.V_long_e));

               v_e = v[2] - (v_l / 2.0);

               do
               { //L_M = U_Ln((R_max+R_min_mag)/R_min_mag);
                  U_rel = (((v[2] - v_e) * (v[2] - v_e) + Vtr * Vtr));
                  Maxw = U_Exp(-1.0 * (v_e * v_e) / (F.V_long_e * F.V_long_e * 2.));
                  norm_m += Maxw;

                  //R_min = F.Z*U_e*U_e/(U_me*(U_rel+F.V_tr_e*F.V_tr_e));
                  if (R_min_mag > R_max)
                     L_M = 0.0;
                  else
                  {
                     L_M = U_Ln(R_max / R_min_mag);
                     //if(R_min > ro_e)L_M = U_Ln((R_max+R_min)/R_min);
                  }

                  d_x_m += Vtr * (Vtr * Vtr - 2.0 * (v[2] - v_e) * (v[2] - v_e)) * L_M * Maxw / (U_rel ^ 2.5);
                  d_s_m += 3.0 * Vtr * Vtr * (v[2] - v_e) * L_M * Maxw / (U_rel * (U_rel ^ 1.5));

                  if (Constant)
                  {
                     if ((U_rel ^ 0.5) > (1.0 * F.V_long_e))
                        d_s_m += 2.0 * (v[2] - v_e) * Maxw / (U_rel ^ 1.5);
                  }
                  v_e -= v_l;
               } while (v_e >= -4.1 * F.V_long_e);

               delta_x_m = d_x_m / norm_m;
               delta_s_m = d_s_m / norm_m;
               //---- 05.06

               if (U_Abs(Vtr) <= 0.01 * F.V_long_e)
               {
                  delta_s_m = 2.0 * L_M * v[2] / (((2.0 * U_pi) ^ 0.5) * F.V_long_e * F.V_long_e * F.V_long_e * U_Exp((v[2] * v[2]) / (2.0 * F.V_long_e * F.V_long_e)));
               }
            }
         }

         delta_x += delta_x_m;
         delta_s += delta_s_m;

         if (Adiabatic)
         { // Adiabatic -------------------------------------------
            if (theta <= F.V_tr_e)
            {
               delta_x_a = Vtr * (2 * (N_col * L_A) / (F.V_tr_e ^ 3));
               delta_s_a = v[2] * (2 * (N_col * L_A) / (F.V_tr_e ^ 2)) / F.V_long_e;

               if (theta >= F.V_long_e)
               {
                  ellipsoid = (Vtr ^ 2) / (F.V_tr_e ^ 2) + ((v[2] ^ 2) / (F.V_long_e ^ 2));
                  if (ellipsoid >= 1.)
                  {
                     t = v[2]() < 0. ? -1. : 1.;
                     delta_s_a = 2 * t * (N_col * L_A) / (F.V_tr_e ^ 2);
                  }
               }
            }
         }

         delta_x += delta_x_a;
         delta_s += delta_s_a;
      }
      Ftr = -delta_x * delta;
      f[2] = -delta_s * delta;
   }
   else
   {
      Ftr = 0;
      f[2] = 0; // total angle is 0
   }
}

void xForce::Parhom(xFrParam F)
{
   delta = 4. * F.n_e * F.Z * F.Z * (U_e ^ 4) / U_me; //Friction constant

   doubleU R_max(m_);  // max impact parameter
   doubleU R_max1(m_); // for max impact parameter (comparisom)
   doubleU R_max2(m_);
   doubleU R_min(m_);                    // min impact parameter
   doubleU ro_e(m_);                     // electron Larmour radius
   doubleU E_e(M_6 * eV_ * 0.5110034);   // rest electron energy
   doubleU theta(m_ / s_);               // total ion velocity
   doubleU theta_e(m_ / s_);             // total electron velocity spread
   doubleU theta_real(m_ / s_);          // total ion velocity & total electron velocity spread
   doubleU L_P;                          // Coloumb logarithm
   doubleU V_3((m_ ^ 3) / (s_ ^ 3));     // cubic of total velocity
   doubleU delta_x((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_s((s_ ^ 2) / (m_ ^ 2)); //

   theta = ((Vtr ^ 2) + (v[2] ^ 2)) ^ 0.5;

   if (theta() > 0.)
   {
      //components of the friction force
      delta_x = 0;
      delta_s = 0;

      theta_e = (F.V_eff_e * F.V_eff_e + F.V_long_e * F.V_long_e) ^ 0.5;
      theta_real = (((theta * theta) + (theta_e * theta_e)) ^ 0.5);

      //maximum and minimum impact parameters constant
      // from Betacool
      // ------------------ 05.06
      R_max = ((theta_real * theta_real * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);

      //         R_max1 = ((3* F.Z /F.n_e)^(1./3.));
      //         if (R_max() < R_max1()) R_max = R_max1; //If Debue sphere does not contain enough electrons R_max is determined by electron density

      R_max2 = theta_real * F.tau;
      if (R_max() > R_max2())
         R_max = R_max2;

      //         R_max = theta_real/( (1/F.tau) + ( (4*U_pi*F.n_e*U_e*U_e/U_me) ^0.5));
      //         R_max = theta /( (1/F.tau) + ((4*U_pi*F.n_e*U_e*U_e/U_me)^0.5));

      //from B		R_min = F.Z*U_e*U_e/(U_me*(theta + F.V_tr_e)*(theta + F.V_tr_e));
      R_min = F.Z * U_e * U_e / (U_me * theta_real * theta_real);

      ro_e = U_me * F.V_tr_e * U_c / (U_e * F.mfield);
      L_P = U_Ln((R_max + R_min + ro_e) / (R_min + ro_e));
      V_3 = ((theta * theta) + (theta_e * theta_e)) * (((theta * theta) + (theta_e * theta_e)) ^ 0.5);

      //            L_P = 2;
      //            V_3 = 1e5;
      delta_x = Vtr * L_P / V_3;
      delta_s = v[2] * L_P / V_3;

      Ftr = -delta_x * delta;
      f[2] = -delta_s * delta;
   }
   else
   {
      Ftr = 0;
      f[2] = 0; // total angle is 0
   }
}

void xForce::Toepffer(xFrParam F)
{
   delta = 4 * U_pi * F.n_e * F.Z * F.Z * (U_e ^ 4) / U_me; //Friction constant
                                                            //   double t;
   doubleU R_max(m_);                                       // max impact parameter
   doubleU R_max1(m_);                                      // for max impact parameter (comparisom)
   doubleU R_max2(m_);
   doubleU R_max_F(m_);                // max impact parameter for fast collisions
   doubleU R_min(m_);                  // min impact parameter
   doubleU R_delta(m_);                // helix parameter
   doubleU ro_e(m_);                   // electron Larmour radius
   doubleU E_e(M_6 * eV_ * 0.5110034); // rest electron energy
   doubleU theta_e(m_ / s_);           // total electron velocity spread
   doubleU theta(m_ / s_);             // total ion velocity
   doubleU L_T;                        // Coloumb logarithm

   // Numerical integration
   doubleU delta_x_f((s_ ^ 2) / (m_ ^ 2)); // fast collisions
   doubleU delta_s_f((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_x_t((s_ ^ 2) / (m_ ^ 2)); // tight helices
   doubleU delta_s_t((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_x_s((s_ ^ 2) / (m_ ^ 2)); // stretched helices
   doubleU delta_s_s((s_ ^ 2) / (m_ ^ 2)); //

   doubleU Maxw;
   doubleU Maxw_t;
   doubleU Maxw_l;

   doubleU d_x((s_) / (m_));           //
   doubleU d_s((s_) / (m_));           //
   doubleU d_x_m((s_ ^ 2) / (m_ ^ 2)); //
   doubleU d_s_m((s_ ^ 2) / (m_ ^ 2)); //

   int i, j, l;

   double dfi; //step over asimuth
   double fi0; //initial asimuth

   doubleU v_l(m_ / s_);   // step over longitudinal electron velocity
   doubleU v_e_l(m_ / s_); // longitudinal electron velocity
   doubleU v_t(m_ / s_);   // step over transverse electron velocity
   doubleU v_e_t(m_ / s_); // transverse electron velocity

   doubleU U_rel((m_ ^ 2) / (s_ ^ 2)); // square of relative velocity

   doubleU norm(m_ / s_); //normalization constant for fast integral
   doubleU norm_m;        //normalization constant for magnetized integrals

   theta = ((Vtr ^ 2) + (v[2] ^ 2)) ^ 0.5;

   if (theta() > 0.)
   {
      //components of the friction force
      delta_x_f = 0;
      delta_s_f = 0;
      delta_x_t = 0;
      delta_s_t = 0;
      delta_x_s = 0;
      delta_s_s = 0;
      //integration steps

      v_t = 3.0 * F.V_tr_e / Tdt;
      v_l = 4.1 * F.V_long_e / Tdl;

      //maximum impact parameter
      R_max = (((theta * theta + F.V_long_e * F.V_long_e) * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);
      R_max1 = ((3 * F.Z / F.n_e) ^ (1. / 3.));
      R_max2 = ((theta * theta + F.V_long_e * F.V_long_e) ^ 0.5) * F.tau;

      //If Debue sphere does not contain enough electrons R_max is determined by electron density
      if (R_max() < R_max1())
         R_max = R_max1;
      if (R_max() > R_max2())
         R_max = R_max2;

      if (TFast)
      {
         R_max_F = R_max;
         dfi = M_PI / Tnfi;
         fi0 = dfi / 2;
         v_t = (F.V_tr_e * 3) / Tdt;
         v_l = (F.V_long_e * 3) / Tdl;

         d_x = 0;
         d_s = 0;
         norm = 0;

         for (i = 1 - Tdl; i < Tdl; i++)
            for (j = 0; j < Tdt; j++)
               for (l = 0; l < Tnfi; l++)
               {
                  U_rel = (((v[2] - (v_l * i)) * (v[2] - (v_l * i)) +
                            (Vtr - ((v_t * j) * cos(dfi * l + fi0))) * (Vtr - ((v_t * j) * cos(dfi * l + fi0))) +
                            (v_t * j) * sin(dfi * l + fi0) * (v_t * j) * sin(dfi * l + fi0)));

                  ro_e = U_me * v_t * j * U_c / (U_e * F.mfield);
                  R_max_F = R_max;
                  if (ro_e() < R_max())
                     R_max_F = ro_e;

                  if (U_rel() != 0.0)
                  {
                     R_min = F.Z * U_e * U_e / (U_me * U_rel);

                     if (R_min < R_max_F)
                     {
                        L_T = 0.5 * U_Ln(1. + ((R_max_F * R_max_F) / (R_min * R_min)));
                     }
                     else
                        L_T = 0.0;

                     Maxw = U_Exp(((v_t * j) * (v_t * j)) / (F.V_tr_e * F.V_tr_e * 2.) +
                                  ((v_l * i) * (v_l * i)) / (F.V_long_e * F.V_long_e * 2.));
                     d_x += (Vtr - ((v_t * j) * cos(dfi * l + fi0))) * L_T * (v_t * j) / ((U_rel ^ 1.5) * Maxw);
                     d_s += (v[2] - (v_l * i)) * L_T * (v_t * j) / ((U_rel ^ 1.5) * Maxw);
                     norm += (v_t * j) / Maxw;
                  }
               }
         delta_x_f = d_x / norm;
         delta_s_f = d_s / norm;
      }

      if (Tight)
      {
         v_t = (F.V_tr_e * 3) / Tdt;
         v_l = (F.V_long_e * 4.1) / N_M;

         d_x_m = 0;
         d_s_m = 0;
         norm_m = 0;
         v_e_l = v[2] + (v_l / 2.0);

         do
         {
            U_rel = (((v[2] - v_e_l) * (v[2] - v_e_l) + Vtr * Vtr));

            Maxw_l = U_Exp(-1.0 * (v_e_l * v_e_l) / (F.V_long_e * F.V_long_e * 2.));
            R_delta = U_me * (U_rel ^ 0.5) * U_c / (U_e * F.mfield);
            L_T = 0;
            // Integral from Coulomb logarithm
            for (j = 0; j < Tdt; j++)
            {
               Maxw_t = U_Exp(-1.0 * ((v_t * j) * (v_t * j)) / (F.V_tr_e * F.V_tr_e * 2.));
               norm_m += Maxw_l * Maxw_t * v_t * j / F.V_tr_e;
               ro_e = U_me * v_t * j * U_c / (U_e * F.mfield);
               R_min = ro_e;
               if (ro_e() < R_delta())
                  R_min = R_delta;
               if (R_max() > R_min())
                  L_T += U_Ln(R_max / R_min) * Maxw_t * v_t * j / F.V_tr_e;
            }
            //
            d_x_m += Vtr * (Vtr * Vtr - (v[2] - v_e_l) * (v[2] - v_e_l)) * L_T * Maxw_l / (U_rel ^ 2.5);
            d_s_m += Vtr * Vtr * (v[2] - v_e_l) * L_T * Maxw_l / (U_rel * (U_rel ^ 1.5));

            v_e_l += v_l;

         } while (v_e_l <= (4.1 * F.V_long_e));

         v_e_l = v[2] - (v_l / 2.0);

         do
         {
            U_rel = (((v[2] - v_e_l) * (v[2] - v_e_l) + Vtr * Vtr));

            Maxw_l = U_Exp(-1.0 * (v_e_l * v_e_l) / (F.V_long_e * F.V_long_e * 2.));
            R_delta = U_me * (U_rel ^ 0.5) * U_c / (U_e * F.mfield);
            L_T = 0;
            // Integral from Coulomb logarithm
            for (j = 0; j < Tdt; j++)
            {
               Maxw_t = U_Exp(-1.0 * ((v_t * j) * (v_t * j)) / (F.V_tr_e * F.V_tr_e * 2.));
               norm_m += Maxw_l * Maxw_t * v_t * j / F.V_tr_e;
               ro_e = U_me * v_t * j * U_c / (U_e * F.mfield);
               R_min = ro_e;
               if (ro_e() < R_delta())
                  R_min = R_delta;
               if (R_max() > R_min())
                  L_T += U_Ln(R_max / R_min) * Maxw_t * v_t * j / F.V_tr_e;
            }
            //
            d_x_m += Vtr * (Vtr * Vtr - (v[2] - v_e_l) * (v[2] - v_e_l)) * L_T * Maxw_l / (U_rel ^ 2.5);
            d_s_m += Vtr * Vtr * (v[2] - v_e_l) * L_T * Maxw_l / (U_rel * (U_rel ^ 1.5));

            v_e_l -= v_l;

         } while (v_e_l >= -4.1 * F.V_long_e);

         delta_x_t = d_x_m / norm_m;
         delta_s_t = d_s_m / norm_m;
      }
      if (Stretched)
      {
         v_t = (F.V_tr_e * 3) / Tdt;
         v_l = (F.V_long_e * 4.1) / N_M;

         d_x_m = 0;
         d_s_m = 0;
         norm_m = 0;
         v_e_l = v[2] + (v_l / 2.0);

         do
         {
            U_rel = (((v[2] - v_e_l) * (v[2] - v_e_l) + Vtr * Vtr));

            Maxw_l = U_Exp(-1.0 * (v_e_l * v_e_l) / (F.V_long_e * F.V_long_e * 2.));
            R_delta = U_me * (U_rel ^ 0.5) * U_c / (U_e * F.mfield);
            L_T = 0;
            // Integral from Coulomb logarithm
            for (j = 0; j < Tdt; j++)
            {
               Maxw_t = U_Exp(-1.0 * ((v_t * j) * (v_t * j)) / (F.V_tr_e * F.V_tr_e * 2.));
               norm_m += Maxw_l * Maxw_t * v_t * j / F.V_tr_e;
               ro_e = U_me * v_t * j * U_c / (U_e * F.mfield);
               R_min = F.Z * U_e * U_e / (U_me * U_rel);
               if (ro_e() > R_min())
                  R_min = ro_e;
               R_max_F = R_delta;
               if (R_max() < R_delta())
                  R_max_F = R_max;
               if (R_min() > R_max())
                  R_min = R_max;
               if (R_max_F() > R_min())
                  L_T += U_Ln(R_max_F / R_min) * Maxw_t * v_t * j / F.V_tr_e;
            }
            //
            d_x_m += Vtr * L_T * Maxw_l / (U_rel ^ 1.5);
            d_s_m += (v[2] - v_e_l) * L_T * Maxw_l / (U_rel ^ 1.5);

            v_e_l += v_l;

         } while (v_e_l <= (4.1 * F.V_long_e));

         v_e_l = v[2] - (v_l / 2.0);

         do
         {
            U_rel = (((v[2] - v_e_l) * (v[2] - v_e_l) + Vtr * Vtr));

            Maxw_l = U_Exp(-1.0 * (v_e_l * v_e_l) / (F.V_long_e * F.V_long_e * 2.));
            R_delta = U_me * (U_rel ^ 0.5) * U_c / (U_e * F.mfield);
            L_T = 0;
            // Integral from Coulomb logarithm
            for (j = 0; j < Tdt; j++)
            {
               Maxw_t = U_Exp(-1.0 * ((v_t * j) * (v_t * j)) / (F.V_tr_e * F.V_tr_e * 2.));
               norm_m += Maxw_l * Maxw_t * v_t * j / F.V_tr_e;
               ro_e = U_me * v_t * j * U_c / (U_e * F.mfield);
               R_min = F.Z * U_e * U_e / (U_me * U_rel);
               if (ro_e() > R_min())
                  R_min = ro_e;
               R_max_F = R_delta;
               if (R_max() < R_delta())
                  R_max_F = R_max;
               if (R_min() > R_max())
                  R_min = R_max;
               if (R_max_F() > R_min())
                  L_T += U_Ln(R_max_F / R_min) * Maxw_t * v_t * j / F.V_tr_e;
            }
            //
            d_x_m += Vtr * L_T * Maxw_l / (U_rel ^ 1.5);
            d_s_m += (v[2] - v_e_l) * L_T * Maxw_l / (U_rel ^ 1.5);

            v_e_l -= v_l;

         } while (v_e_l >= -4.1 * F.V_long_e);

         delta_x_s = d_x_m / norm_m;
         delta_s_s = d_s_m / norm_m;
      }

      Ftr = -(delta_x_f + delta_x_t + delta_x_s) * delta;
      f[2] = -(delta_s_f + delta_s_t + delta_s_s) * delta;
   }
   else
   {
      Ftr = 0;
      f[2] = 0; // total angle is 0
   }
}

void xForce::D4(xFrParam F)
{
   delta = 4 * U_pi * F.n_e * F.Z * F.Z * (U_e ^ 4) / U_me; //Friction constant

   doubleU R_max(m_);  // max impact parameter
   doubleU R_max1(m_); // for max impact parameter (comparisom)
   doubleU R_max2(m_);
   doubleU R_min(m_); // min impact parameter

   doubleU theta(m_ / s_);   // total ion velocity
   doubleU theta_e(m_ / s_); //total electron velocity

   doubleU norm(U1_);

   doubleU L_C; // Coloumb logarithms

   doubleU Maxw;

   doubleU delta_x((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_y((s_ ^ 2) / (m_ ^ 2)); //
   doubleU delta_s((s_ ^ 2) / (m_ ^ 2)); //

   int i, j, l;

   doubleU v_l(m_ / s_); // step over longitudinal electron velocity
   doubleU v_x(m_ / s_); // step over horizontal electron velocity
   doubleU v_y(m_ / s_); // step over vertical electron velocity

   doubleU U_rel((m_ ^ 2) / (s_ ^ 2)); // square of relative velocity

   theta = ((v[0] ^ 2) + (v[1] ^ 2) + (v[2] ^ 2)) ^ 0.5;
   theta_e = ((F.V_tr_x ^ 2) + (F.V_tr_y ^ 2) + (F.V_long_e ^ 2)) ^ 0.5;

   if (theta() > 0.)
   {
      //steps over electron velocities
      v_x = (F.V_tr_x * 3) / D3dx;
      v_y = (F.V_tr_y * 3) / D3dy;
      v_l = (F.V_long_e * 3) / D3dl;
      //components of the friction force
      delta_x = 0;
      delta_y = 0;
      delta_s = 0;
      //intgral over electron distribution (normalization factor)
      norm = 0;

      //maximum impact parameter constant
      R_max1 = (3 * F.Z / F.n_e) ^ (1. / 3.);
      R_max2 = theta * F.tau;

      if (theta > theta_e)
         R_max = (theta * theta * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5;
      else
         R_max = (theta_e * theta_e * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5;

      if (R_max < R_max1)
         R_max = R_max1; // choose max
      if (R_max > R_max2)
         R_max = R_max2; // choose min

      for (i = 1 - D3dl; i < D3dl; i++)
         for (j = 1 - D3dx; j < D3dx; j++)
            for (l = 1 - D3dy; l < D3dy; l++)
            {
               U_rel = ((v[2] - (v_l * i)) * (v[2] - (v_l * i))) +
                       ((v[0] - (v_x * j)) * (v[0] - (v_x * j))) +
                       ((v[1] - (v_y * l)) * (v[1] - (v_y * l)));

               Maxw = U_Exp(((v_x * j) * (v_x * j)) / (F.V_tr_x * F.V_tr_x * 2.) +
                            ((v_y * l) * (v_y * l)) / (F.V_tr_y * F.V_tr_y * 2.) +
                            ((v_l * i) * (v_l * i)) / (F.V_long_e * F.V_long_e * 2.));

               if (U_rel() != 0.0)
               {
                  // minimum impact parameter
                  R_min = F.Z * U_e * U_e / (U_me * U_rel);
                  if (F.undulator && R_min < F.r_0)
                     R_min = F.r_0;

                  if (R_max() > R_min())
                  {
                     L_C = U_Ln(R_max / R_min);
                     delta_x += (v[0] - (v_x * j)) * L_C / ((U_rel ^ 1.5) * Maxw);
                     delta_y += (v[1] - (v_y * l)) * L_C / ((U_rel ^ 1.5) * Maxw);
                     delta_s += (v[2] - (v_l * i)) * L_C / ((U_rel ^ 1.5) * Maxw);
                  }
               }
               norm += 1.0 / Maxw;
            }
      delta_x /= norm;
      delta_y /= norm;
      delta_s /= norm;

      f[0] = -delta_x * delta;
      f[1] = -delta_y * delta;
      f[2] = -delta_s * delta;
   }
   else
   {
      f[0] = 0;
      f[1] = 0;
      f[2] = 0; // total angle is 0
   }
}

void xForce::Linear(xFrParam F)
{
   delta = 2 * U_pi * F.n_e * F.Z * F.Z * (U_e ^ 4) / U_me; //Friction constant
   doubleU R_max(m_);                                       // max impact parameter
   doubleU R_max1(m_);                                      // for max impact parameter (comparisom)
   doubleU R_min(m_);                                       // min impact parameter
   doubleU R_f(m_);                                         // intermediate impact parameter
   doubleU ro_e(m_);                                        // electron Larmour radius
   doubleU E_e(M_6 * eV_ * 0.5110034);                      // rest electron energy
   doubleU theta_e(m_ / s_);                                // total electron velocity spread
   doubleU theta(m_ / s_);                                  // total ion velocity
   doubleU L_M, L_A, L_F;                                   // Coloumb logarithms
   doubleU N_col;                                           // number of collisions
   doubleU delta_x((s_ ^ 2) / (m_ ^ 2));                    //
   doubleU delta_s((s_ ^ 2) / (m_ ^ 2));                    //
   doubleU ellipsoid;                                       // constant determined region of velocity

   //Components of the ion velocity
   theta = ((Vtr ^ 2) + (v[2] ^ 2)) ^ 0.5;

   if (theta() > 0.)
   {
      //components of the friction force
      delta_x = 0;
      delta_s = 0;
      //maximum and minimum impact parameters constant
      R_max = ((theta * theta * U_me / (4 * U_pi * F.n_e * U_e * U_e)) ^ 0.5);
      R_max1 = ((3 * F.Z / F.n_e) ^ (1. / 3.));
      if (R_max() > R_max1())
         R_max = R_max1;
      R_min = F.Z * U_e * U_e / (U_me * (theta + F.V_tr_e) * (theta + F.V_tr_e));

      theta_e = ((F.V_tr_e * F.V_tr_e + F.V_long_e * F.V_long_e) ^ 0.5);
      ro_e = U_me * F.V_tr_e * U_c / (U_e * F.mfield);
      R_f = ro_e * (theta + F.V_long_e) / F.V_tr_e;
      L_M = 0;
      L_A = 0;
      L_F = 0;
      if (R_max > (2 * ro_e))
         L_M = U_Ln(R_max / (2 * ro_e));
      if ((2 * ro_e) > R_f)
         L_A = U_Ln(2 * ro_e / R_f);
      if (R_f > R_min)
         L_F = U_Ln(R_f / R_min);

      N_col = 1 + F.V_tr_e / (U_pi * (theta + F.V_long_e));

      delta_x = Vtr * (2 * (L_F + N_col * L_A) / (F.V_tr_e ^ 3) + (L_M / (F.V_long_e ^ 3)));
      delta_s = v[2] * (2 * (L_F + N_col * L_A) / (F.V_tr_e ^ 2) + (L_M / (F.V_long_e ^ 2))) / F.V_long_e;

      Ftr = -delta_x * delta;
      f[2] = -delta_s * delta;
   }
   else
   {
      Ftr = 0;
      f[2] = 0; // total angle is 0
   }
}
//----------------------------------------------------------------------------
//Error   integral from R.Zeinalov
//--------------------------------------------------------------------------
double xForce::phi_f(double x)
{
   return exp(-x * x / 2);
}

double xForce::Erf(double x)
{

   int spl;

   if (x <= 1)
      spl = 200;
   else if (x > 1 && x < 3)
      spl = 200;
   else if (x >= 3 && x < 6)
      spl = 400;
   else if (x >= 6 && x < 12)
      spl = 800;
   else if (x >= 12 && x < 24)
      spl = 1600;
   else if (x >= 24 && x < 48)
      spl = 3200;
   else if (x >= 48 && x < 96)
      spl = 6400;
   else
      spl = 12800;

   double step = x;
   double cx = 0;
   double integ = 0;
   double linteg = 0;
   double tmp = 0;

   if (x == 0)
   {
      return 0;
   }

   do
   {
      linteg = integ;
      integ = 0;
      spl *= 2;
      step = x / spl;

      cx = 0;
      while (cx <= x - step)
      {
         integ += phi_f(cx);
         cx += step;
      }
      integ *= 2;

      tmp = 0;
      cx = 3 * step / 2;
      while (cx <= x - 0.5 * step)
      {
         tmp += phi_f(cx);
         cx += step;
      }
      tmp *= 4;
      integ += tmp;

      integ += phi_f(0) + phi_f(x);
      integ *= step / 6;
   } while (linteg == 0);

   return integ;
}

double xForce::phi(double x)
{
   double integ;
   double pi = 3.14159265;

   integ = sqrt(2 / pi) * (Erf(x) - x * phi_f(x));

   return integ;
}
