#include "NSEDualGradient3D.h"
/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_0::UFC_NSEDualGradient3DLinearForm_finite_element_0_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_0::~UFC_NSEDualGradient3DLinearForm_finite_element_0_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_0::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_0::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_0::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_0::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_0::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_0::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_0();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_1::UFC_NSEDualGradient3DLinearForm_finite_element_0_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_1::~UFC_NSEDualGradient3DLinearForm_finite_element_0_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_1::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_1::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_1::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_1::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_1::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_1();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_2::UFC_NSEDualGradient3DLinearForm_finite_element_0_2() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_2::~UFC_NSEDualGradient3DLinearForm_finite_element_0_2()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_2::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_2::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_2::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_2::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_2::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_2::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_2::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_2::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_2::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_2::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_2::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_2::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_2::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_2::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_2();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_3::UFC_NSEDualGradient3DLinearForm_finite_element_0_3() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_3::~UFC_NSEDualGradient3DLinearForm_finite_element_0_3()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_3::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_3::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_3::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_3::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_3::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_3::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_3::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_3::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_3::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_3::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_3::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_3::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_3::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_3::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_3();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_4::UFC_NSEDualGradient3DLinearForm_finite_element_0_4() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_4::~UFC_NSEDualGradient3DLinearForm_finite_element_0_4()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_4::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_4::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_4::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_4::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_4::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_4::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_4::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_4::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_4::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_4::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_4::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_4::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_4::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_4::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_4();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_5::UFC_NSEDualGradient3DLinearForm_finite_element_0_5() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_5::~UFC_NSEDualGradient3DLinearForm_finite_element_0_5()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_5::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_5::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_5::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_5::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_5::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_5::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_5::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_5::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_5::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_5::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_5::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_5::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_5::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_5::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_5();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_6::UFC_NSEDualGradient3DLinearForm_finite_element_0_6() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_6::~UFC_NSEDualGradient3DLinearForm_finite_element_0_6()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_6::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_6::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_6::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_6::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_6::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_6::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_6::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_6::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_6::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_6::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_6::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_6::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_6::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_6::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_6();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_7::UFC_NSEDualGradient3DLinearForm_finite_element_0_7() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_7::~UFC_NSEDualGradient3DLinearForm_finite_element_0_7()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_7::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_7::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_7::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_7::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_7::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_7::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_7::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_7::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_7::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_7::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_7::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_7::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_7::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_7::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_7();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_8::UFC_NSEDualGradient3DLinearForm_finite_element_0_8() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0_8::~UFC_NSEDualGradient3DLinearForm_finite_element_0_8()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0_8::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0_8::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_8::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_8::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_8::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_8::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_8::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_8::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0_8::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0_8::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0_8::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0_8::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0_8::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0_8::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_0_8();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_0::UFC_NSEDualGradient3DLinearForm_finite_element_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_0::~UFC_NSEDualGradient3DLinearForm_finite_element_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_0::signature() const
{
    return "Mixed finite element: [Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron]";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0::space_dimension() const
{
    return 9;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0::value_dimension(unsigned int i) const
{
    return 9;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    values[4] = 0;
    values[5] = 0;
    values[6] = 0;
    values[7] = 0;
    values[8] = 0;
    
    if (0 <= i && i <= 0)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[0] = coeff0_0*basisvalue0;
    }
    
    if (1 <= i && i <= 1)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 1;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[1] = coeff0_0*basisvalue0;
    }
    
    if (2 <= i && i <= 2)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 2;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[2] = coeff0_0*basisvalue0;
    }
    
    if (3 <= i && i <= 3)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 3;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[3] = coeff0_0*basisvalue0;
    }
    
    if (4 <= i && i <= 4)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[4] = coeff0_0*basisvalue0;
    }
    
    if (5 <= i && i <= 5)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 5;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[5] = coeff0_0*basisvalue0;
    }
    
    if (6 <= i && i <= 6)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 6;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[6] = coeff0_0*basisvalue0;
    }
    
    if (7 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 7;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[7] = coeff0_0*basisvalue0;
    }
    
    if (8 <= i && i <= 8)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 8;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[8] = coeff0_0*basisvalue0;
    }
    
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 9*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 0)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (1 <= i && i <= 1)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 1;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (2 <= i && i <= 2)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 2;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[2*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (3 <= i && i <= 3)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 3;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[3*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (4 <= i && i <= 4)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[4*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (5 <= i && i <= 5)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 5;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[5*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (6 <= i && i <= 6)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 6;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[6*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (7 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 7;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[7*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (8 <= i && i <= 8)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 8;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_z_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
      // Table(s) of coefficients
      const static double coefficients0[1][1] =   \
      {{1.15470053837925}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[1][1] =   \
      {{0}};
    
      const static double dmats1[1][1] =   \
      {{0}};
    
      const static double dmats2[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[8*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[9][1][3] = {{{0.25, 0.25, 0.25}}, {{0.25, 0.25, 0.25}}, {{0.25, 0.25, 0.25}}, {{0.25, 0.25, 0.25}}, {{0.25, 0.25, 0.25}}, {{0.25, 0.25, 0.25}}, {{0.25, 0.25, 0.25}}, {{0.25, 0.25, 0.25}}, {{0.25, 0.25, 0.25}}};
    const static double W[9][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[9][1][9] = {{{1, 0, 0, 0, 0, 0, 0, 0, 0}}, {{0, 1, 0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 1, 0, 0, 0, 0, 0, 0}}, {{0, 0, 0, 1, 0, 0, 0, 0, 0}}, {{0, 0, 0, 0, 1, 0, 0, 0, 0}}, {{0, 0, 0, 0, 0, 1, 0, 0, 0}}, {{0, 0, 0, 0, 0, 0, 1, 0, 0}}, {{0, 0, 0, 0, 0, 0, 0, 1, 0}}, {{0, 0, 0, 0, 0, 0, 0, 0, 1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[9];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 9; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_0::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[9] = dof_values[0];
    vertex_values[18] = dof_values[0];
    vertex_values[27] = dof_values[0];
    // Evaluate at vertices and use affine mapping
    vertex_values[1] = dof_values[1];
    vertex_values[10] = dof_values[1];
    vertex_values[19] = dof_values[1];
    vertex_values[28] = dof_values[1];
    // Evaluate at vertices and use affine mapping
    vertex_values[2] = dof_values[2];
    vertex_values[11] = dof_values[2];
    vertex_values[20] = dof_values[2];
    vertex_values[29] = dof_values[2];
    // Evaluate at vertices and use affine mapping
    vertex_values[3] = dof_values[3];
    vertex_values[12] = dof_values[3];
    vertex_values[21] = dof_values[3];
    vertex_values[30] = dof_values[3];
    // Evaluate at vertices and use affine mapping
    vertex_values[4] = dof_values[4];
    vertex_values[13] = dof_values[4];
    vertex_values[22] = dof_values[4];
    vertex_values[31] = dof_values[4];
    // Evaluate at vertices and use affine mapping
    vertex_values[5] = dof_values[5];
    vertex_values[14] = dof_values[5];
    vertex_values[23] = dof_values[5];
    vertex_values[32] = dof_values[5];
    // Evaluate at vertices and use affine mapping
    vertex_values[6] = dof_values[6];
    vertex_values[15] = dof_values[6];
    vertex_values[24] = dof_values[6];
    vertex_values[33] = dof_values[6];
    // Evaluate at vertices and use affine mapping
    vertex_values[7] = dof_values[7];
    vertex_values[16] = dof_values[7];
    vertex_values[25] = dof_values[7];
    vertex_values[34] = dof_values[7];
    // Evaluate at vertices and use affine mapping
    vertex_values[8] = dof_values[8];
    vertex_values[17] = dof_values[8];
    vertex_values[26] = dof_values[8];
    vertex_values[35] = dof_values[8];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_0::num_sub_elements() const
{
    return 9;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_0::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_0();
      break;
    case 1:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_1();
      break;
    case 2:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_2();
      break;
    case 3:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_3();
      break;
    case 4:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_4();
      break;
    case 5:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_5();
      break;
    case 6:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_6();
      break;
    case 7:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_7();
      break;
    case 8:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0_8();
      break;
    }
    return 0;
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_1_0::UFC_NSEDualGradient3DLinearForm_finite_element_1_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_1_0::~UFC_NSEDualGradient3DLinearForm_finite_element_1_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_1_0::signature() const
{
    return "Lagrange finite element of degree 1 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_1_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_0::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_0::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_10_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
    // Table(s) of coefficients
    const static double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0, 0, 0.223606797749979}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_0::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_10_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
    // Table(s) of coefficients
    const static double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0, 0, 0.223606797749979}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[4][4] = \
    {{0, 0, 0, 0},
    {6.32455532033676, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0}};
    
    const static double dmats1[4][4] = \
    {{0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0},
    {5.47722557505166, 0, 0, 0},
    {0, 0, 0, 0}};
    
    const static double dmats2[4][4] = \
    {{0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0},
    {1.82574185835055, 0, 0, 0},
    {5.16397779494322, 0, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff0_3 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff0_3 = coefficients0[dof][3];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff0_3 = new_coeff0_3;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0] + coeff0_1*dmats2[1][0] + coeff0_2*dmats2[2][0] + coeff0_3*dmats2[3][0];
          new_coeff0_1 = coeff0_0*dmats2[0][1] + coeff0_1*dmats2[1][1] + coeff0_2*dmats2[2][1] + coeff0_3*dmats2[3][1];
          new_coeff0_2 = coeff0_0*dmats2[0][2] + coeff0_1*dmats2[1][2] + coeff0_2*dmats2[2][2] + coeff0_3*dmats2[3][2];
          new_coeff0_3 = coeff0_0*dmats2[0][3] + coeff0_1*dmats2[1][3] + coeff0_2*dmats2[2][3] + coeff0_3*dmats2[3][3];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_1_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[4][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[4][1] = {{1}, {1}, {1}, {1}};
    const static double D[4][1][1] = {{{1}}, {{1}}, {{1}}, {{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_1_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_1_0::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
    vertex_values[3] = dof_values[3];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_1_0::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_1_0();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_1_1::UFC_NSEDualGradient3DLinearForm_finite_element_1_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_1_1::~UFC_NSEDualGradient3DLinearForm_finite_element_1_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_1_1::signature() const
{
    return "Lagrange finite element of degree 1 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_1_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_1::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_1::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_10_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
    // Table(s) of coefficients
    const static double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0, 0, 0.223606797749979}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_1::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_10_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
    // Table(s) of coefficients
    const static double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0, 0, 0.223606797749979}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[4][4] = \
    {{0, 0, 0, 0},
    {6.32455532033676, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0}};
    
    const static double dmats1[4][4] = \
    {{0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0},
    {5.47722557505166, 0, 0, 0},
    {0, 0, 0, 0}};
    
    const static double dmats2[4][4] = \
    {{0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0},
    {1.82574185835055, 0, 0, 0},
    {5.16397779494322, 0, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff0_3 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff0_3 = coefficients0[dof][3];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff0_3 = new_coeff0_3;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0] + coeff0_1*dmats2[1][0] + coeff0_2*dmats2[2][0] + coeff0_3*dmats2[3][0];
          new_coeff0_1 = coeff0_0*dmats2[0][1] + coeff0_1*dmats2[1][1] + coeff0_2*dmats2[2][1] + coeff0_3*dmats2[3][1];
          new_coeff0_2 = coeff0_0*dmats2[0][2] + coeff0_1*dmats2[1][2] + coeff0_2*dmats2[2][2] + coeff0_3*dmats2[3][2];
          new_coeff0_3 = coeff0_0*dmats2[0][3] + coeff0_1*dmats2[1][3] + coeff0_2*dmats2[2][3] + coeff0_3*dmats2[3][3];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_1_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[4][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[4][1] = {{1}, {1}, {1}, {1}};
    const static double D[4][1][1] = {{{1}}, {{1}}, {{1}}, {{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_1_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_1_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
    vertex_values[3] = dof_values[3];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_1_1::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_1_1();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_1_2::UFC_NSEDualGradient3DLinearForm_finite_element_1_2() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_1_2::~UFC_NSEDualGradient3DLinearForm_finite_element_1_2()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_1_2::signature() const
{
    return "Lagrange finite element of degree 1 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_1_2::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_2::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_2::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_2::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_2::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_10_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
    // Table(s) of coefficients
    const static double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0, 0, 0.223606797749979}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_2::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_2::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_10_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
    // Table(s) of coefficients
    const static double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0, 0, 0.223606797749979}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[4][4] = \
    {{0, 0, 0, 0},
    {6.32455532033676, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0}};
    
    const static double dmats1[4][4] = \
    {{0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0},
    {5.47722557505166, 0, 0, 0},
    {0, 0, 0, 0}};
    
    const static double dmats2[4][4] = \
    {{0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0},
    {1.82574185835055, 0, 0, 0},
    {5.16397779494322, 0, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff0_3 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff0_3 = coefficients0[dof][3];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff0_3 = new_coeff0_3;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0] + coeff0_1*dmats2[1][0] + coeff0_2*dmats2[2][0] + coeff0_3*dmats2[3][0];
          new_coeff0_1 = coeff0_0*dmats2[0][1] + coeff0_1*dmats2[1][1] + coeff0_2*dmats2[2][1] + coeff0_3*dmats2[3][1];
          new_coeff0_2 = coeff0_0*dmats2[0][2] + coeff0_1*dmats2[1][2] + coeff0_2*dmats2[2][2] + coeff0_3*dmats2[3][2];
          new_coeff0_3 = coeff0_0*dmats2[0][3] + coeff0_1*dmats2[1][3] + coeff0_2*dmats2[2][3] + coeff0_3*dmats2[3][3];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1_2::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_1_2::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[4][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[4][1] = {{1}, {1}, {1}, {1}};
    const static double D[4][1][1] = {{{1}}, {{1}}, {{1}}, {{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_1_2::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_1_2::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
    vertex_values[3] = dof_values[3];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1_2::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_1_2::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_1_2();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_1::UFC_NSEDualGradient3DLinearForm_finite_element_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_1::~UFC_NSEDualGradient3DLinearForm_finite_element_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_1::signature() const
{
    return "Mixed finite element: [Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron]";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1::space_dimension() const
{
    return 12;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1::value_dimension(unsigned int i) const
{
    return 3;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    
    if (0 <= i && i <= 3)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_z_0 = 1;
      const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_1_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
      const double psitilde_cs_00_1 = 2*z + 1;
      const double psitilde_cs_01_0 = 1;
      const double psitilde_cs_10_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
      const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
      const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
      const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
      // Table(s) of coefficients
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[0] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_z_0 = 1;
      const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_1_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
      const double psitilde_cs_00_1 = 2*z + 1;
      const double psitilde_cs_01_0 = 1;
      const double psitilde_cs_10_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
      const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
      const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
      const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
      // Table(s) of coefficients
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[1] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
    if (8 <= i && i <= 11)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 8;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_z_0 = 1;
      const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_1_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
      const double psitilde_cs_00_1 = 2*z + 1;
      const double psitilde_cs_01_0 = 1;
      const double psitilde_cs_10_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
      const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
      const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
      const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
      // Table(s) of coefficients
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[2] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 3*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 3)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_z_0 = 1;
      const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_1_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
      const double psitilde_cs_00_1 = 2*z + 1;
      const double psitilde_cs_01_0 = 1;
      const double psitilde_cs_10_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
      const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
      const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
      const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
      // Table(s) of coefficients
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {1.82574185835055, 0, 0, 0},
      {5.16397779494322, 0, 0, 0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
      double coeff0_1 = 0;
      double coeff0_2 = 0;
      double coeff0_3 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
      double new_coeff0_1 = 0;
      double new_coeff0_2 = 0;
      double new_coeff0_3 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
        new_coeff0_1 = coefficients0[dof][1];
        new_coeff0_2 = coefficients0[dof][2];
        new_coeff0_3 = coefficients0[dof][3];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
          coeff0_1 = new_coeff0_1;
          coeff0_2 = new_coeff0_2;
          coeff0_3 = new_coeff0_3;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0];
            new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1];
            new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2];
            new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0];
            new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1];
            new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2];
            new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0] + coeff0_1*dmats2[1][0] + coeff0_2*dmats2[2][0] + coeff0_3*dmats2[3][0];
            new_coeff0_1 = coeff0_0*dmats2[0][1] + coeff0_1*dmats2[1][1] + coeff0_2*dmats2[2][1] + coeff0_3*dmats2[3][1];
            new_coeff0_2 = coeff0_0*dmats2[0][2] + coeff0_1*dmats2[1][2] + coeff0_2*dmats2[2][2] + coeff0_3*dmats2[3][2];
            new_coeff0_3 = coeff0_0*dmats2[0][3] + coeff0_1*dmats2[1][3] + coeff0_2*dmats2[2][3] + coeff0_3*dmats2[3][3];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_z_0 = 1;
      const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_1_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
      const double psitilde_cs_00_1 = 2*z + 1;
      const double psitilde_cs_01_0 = 1;
      const double psitilde_cs_10_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
      const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
      const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
      const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
      // Table(s) of coefficients
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {1.82574185835055, 0, 0, 0},
      {5.16397779494322, 0, 0, 0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
      double coeff0_1 = 0;
      double coeff0_2 = 0;
      double coeff0_3 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
      double new_coeff0_1 = 0;
      double new_coeff0_2 = 0;
      double new_coeff0_3 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
        new_coeff0_1 = coefficients0[dof][1];
        new_coeff0_2 = coefficients0[dof][2];
        new_coeff0_3 = coefficients0[dof][3];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
          coeff0_1 = new_coeff0_1;
          coeff0_2 = new_coeff0_2;
          coeff0_3 = new_coeff0_3;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0];
            new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1];
            new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2];
            new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0];
            new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1];
            new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2];
            new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0] + coeff0_1*dmats2[1][0] + coeff0_2*dmats2[2][0] + coeff0_3*dmats2[3][0];
            new_coeff0_1 = coeff0_0*dmats2[0][1] + coeff0_1*dmats2[1][1] + coeff0_2*dmats2[2][1] + coeff0_3*dmats2[3][1];
            new_coeff0_2 = coeff0_0*dmats2[0][2] + coeff0_1*dmats2[1][2] + coeff0_2*dmats2[2][2] + coeff0_3*dmats2[3][2];
            new_coeff0_3 = coeff0_0*dmats2[0][3] + coeff0_1*dmats2[1][3] + coeff0_2*dmats2[2][3] + coeff0_3*dmats2[3][3];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (8 <= i && i <= 11)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 8;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_z_0 = 1;
      const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_1_0 = 1;
    
      // Compute psitilde_cs
      const double psitilde_cs_00_0 = 1;
      const double psitilde_cs_00_1 = 2*z + 1;
      const double psitilde_cs_01_0 = 1;
      const double psitilde_cs_10_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
      const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
      const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
      const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
      // Table(s) of coefficients
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {1.82574185835055, 0, 0, 0},
      {5.16397779494322, 0, 0, 0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
      double coeff0_1 = 0;
      double coeff0_2 = 0;
      double coeff0_3 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
      double new_coeff0_1 = 0;
      double new_coeff0_2 = 0;
      double new_coeff0_3 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
        new_coeff0_1 = coefficients0[dof][1];
        new_coeff0_2 = coefficients0[dof][2];
        new_coeff0_3 = coefficients0[dof][3];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
          coeff0_1 = new_coeff0_1;
          coeff0_2 = new_coeff0_2;
          coeff0_3 = new_coeff0_3;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0];
            new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1];
            new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2];
            new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0];
            new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1];
            new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2];
            new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3];
          }
          if(combinations[deriv_num][j] == 2)
          {
            new_coeff0_0 = coeff0_0*dmats2[0][0] + coeff0_1*dmats2[1][0] + coeff0_2*dmats2[2][0] + coeff0_3*dmats2[3][0];
            new_coeff0_1 = coeff0_0*dmats2[0][1] + coeff0_1*dmats2[1][1] + coeff0_2*dmats2[2][1] + coeff0_3*dmats2[3][1];
            new_coeff0_2 = coeff0_0*dmats2[0][2] + coeff0_1*dmats2[1][2] + coeff0_2*dmats2[2][2] + coeff0_3*dmats2[3][2];
            new_coeff0_3 = coeff0_0*dmats2[0][3] + coeff0_1*dmats2[1][3] + coeff0_2*dmats2[2][3] + coeff0_3*dmats2[3][3];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[2*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[12][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[12][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[12][1][3] = {{{1, 0, 0}}, {{1, 0, 0}}, {{1, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 1, 0}}, {{0, 1, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 1}}, {{0, 0, 1}}, {{0, 0, 1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[3];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 3; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[3] = dof_values[1];
    vertex_values[6] = dof_values[2];
    vertex_values[9] = dof_values[3];
    // Evaluate at vertices and use affine mapping
    vertex_values[1] = dof_values[4];
    vertex_values[4] = dof_values[5];
    vertex_values[7] = dof_values[6];
    vertex_values[10] = dof_values[7];
    // Evaluate at vertices and use affine mapping
    vertex_values[2] = dof_values[8];
    vertex_values[5] = dof_values[9];
    vertex_values[8] = dof_values[10];
    vertex_values[11] = dof_values[11];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_1::num_sub_elements() const
{
    return 3;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_1::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_1_0();
      break;
    case 1:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_1_1();
      break;
    case 2:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_1_2();
      break;
    }
    return 0;
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_finite_element_2::UFC_NSEDualGradient3DLinearForm_finite_element_2() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_finite_element_2::~UFC_NSEDualGradient3DLinearForm_finite_element_2()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* UFC_NSEDualGradient3DLinearForm_finite_element_2::signature() const
{
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return the cell shape
ufc::shape UFC_NSEDualGradient3DLinearForm_finite_element_2::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_2::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_2::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_2::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_2::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_2::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_2::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void UFC_NSEDualGradient3DLinearForm_finite_element_2::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double UFC_NSEDualGradient3DLinearForm_finite_element_2::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void UFC_NSEDualGradient3DLinearForm_finite_element_2::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void UFC_NSEDualGradient3DLinearForm_finite_element_2::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_finite_element_2::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* UFC_NSEDualGradient3DLinearForm_finite_element_2::create_sub_element(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_finite_element_2();
}

/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_0::UFC_NSEDualGradient3DLinearForm_dof_map_0_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_0::~UFC_NSEDualGradient3DLinearForm_dof_map_0_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_0::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_0::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_0::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_0::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_0::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_0();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_1::UFC_NSEDualGradient3DLinearForm_dof_map_0_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_1::~UFC_NSEDualGradient3DLinearForm_dof_map_0_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_1::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_1::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_1::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_1::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_1::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_1::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_1();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_2::UFC_NSEDualGradient3DLinearForm_dof_map_0_2() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_2::~UFC_NSEDualGradient3DLinearForm_dof_map_0_2()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_2::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_2::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_2::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_2::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_2::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_2::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_2::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_2::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_2::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_2::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_2::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_2::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_2::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_2::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_2::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_2::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_2();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_3::UFC_NSEDualGradient3DLinearForm_dof_map_0_3() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_3::~UFC_NSEDualGradient3DLinearForm_dof_map_0_3()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_3::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_3::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_3::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_3::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_3::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_3::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_3::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_3::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_3::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_3::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_3::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_3::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_3::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_3::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_3::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_3::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_3();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_4::UFC_NSEDualGradient3DLinearForm_dof_map_0_4() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_4::~UFC_NSEDualGradient3DLinearForm_dof_map_0_4()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_4::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_4::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_4::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_4::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_4::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_4::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_4::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_4::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_4::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_4::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_4::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_4::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_4::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_4::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_4::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_4::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_4();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_5::UFC_NSEDualGradient3DLinearForm_dof_map_0_5() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_5::~UFC_NSEDualGradient3DLinearForm_dof_map_0_5()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_5::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_5::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_5::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_5::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_5::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_5::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_5::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_5::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_5::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_5::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_5::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_5::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_5::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_5::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_5::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_5::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_5();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_6::UFC_NSEDualGradient3DLinearForm_dof_map_0_6() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_6::~UFC_NSEDualGradient3DLinearForm_dof_map_0_6()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_6::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_6::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_6::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_6::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_6::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_6::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_6::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_6::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_6::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_6::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_6::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_6::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_6::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_6::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_6::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_6::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_6();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_7::UFC_NSEDualGradient3DLinearForm_dof_map_0_7() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_7::~UFC_NSEDualGradient3DLinearForm_dof_map_0_7()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_7::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_7::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_7::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_7::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_7::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_7::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_7::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_7::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_7::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_7::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_7::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_7::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_7::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_7::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_7::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_7::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_7();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_8::UFC_NSEDualGradient3DLinearForm_dof_map_0_8() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0_8::~UFC_NSEDualGradient3DLinearForm_dof_map_0_8()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0_8::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_8::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0_8::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_8::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0_8::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_8::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_8::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_8::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_8::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_8::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_8::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0_8::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0_8::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0_8::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0_8::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0_8::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_0_8();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_0::UFC_NSEDualGradient3DLinearForm_dof_map_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_0::~UFC_NSEDualGradient3DLinearForm_dof_map_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_0::signature() const
{
    return "FFC dof map for Mixed finite element: [Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron, Discontinuous Lagrange finite element of degree 0 on a tetrahedron]";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_0::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 9*m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0::local_dimension() const
{
    return 9;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
    unsigned int offset = m.num_entities[3];
    dofs[1] = offset + c.entity_indices[3][0];
    offset = offset + m.num_entities[3];
    dofs[2] = offset + c.entity_indices[3][0];
    offset = offset + m.num_entities[3];
    dofs[3] = offset + c.entity_indices[3][0];
    offset = offset + m.num_entities[3];
    dofs[4] = offset + c.entity_indices[3][0];
    offset = offset + m.num_entities[3];
    dofs[5] = offset + c.entity_indices[3][0];
    offset = offset + m.num_entities[3];
    dofs[6] = offset + c.entity_indices[3][0];
    offset = offset + m.num_entities[3];
    dofs[7] = offset + c.entity_indices[3][0];
    offset = offset + m.num_entities[3];
    dofs[8] = offset + c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
    coordinates[1][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[1][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[1][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
    coordinates[2][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[2][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[2][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
    coordinates[3][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[3][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[3][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
    coordinates[4][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[4][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[4][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
    coordinates[5][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[5][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[5][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
    coordinates[6][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[6][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[6][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
    coordinates[7][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[7][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[7][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
    coordinates[8][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[8][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[8][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_0::num_sub_dof_maps() const
{
    return 9;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_0::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_0();
      break;
    case 1:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_1();
      break;
    case 2:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_2();
      break;
    case 3:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_3();
      break;
    case 4:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_4();
      break;
    case 5:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_5();
      break;
    case 6:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_6();
      break;
    case 7:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_7();
      break;
    case 8:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0_8();
      break;
    }
    return 0;
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_1_0::UFC_NSEDualGradient3DLinearForm_dof_map_1_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_1_0::~UFC_NSEDualGradient3DLinearForm_dof_map_1_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_1_0::signature() const
{
    return "FFC dof map for Lagrange finite element of degree 1 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_1_0::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return true;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return false;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_1_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_1_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_0::local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_1_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_1_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = x[0][0];
    coordinates[0][1] = x[0][1];
    coordinates[0][2] = x[0][2];
    coordinates[1][0] = x[1][0];
    coordinates[1][1] = x[1][1];
    coordinates[1][2] = x[1][2];
    coordinates[2][0] = x[2][0];
    coordinates[2][1] = x[2][1];
    coordinates[2][2] = x[2][2];
    coordinates[3][0] = x[3][0];
    coordinates[3][1] = x[3][1];
    coordinates[3][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_1_0::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_1_0();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_1_1::UFC_NSEDualGradient3DLinearForm_dof_map_1_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_1_1::~UFC_NSEDualGradient3DLinearForm_dof_map_1_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_1_1::signature() const
{
    return "FFC dof map for Lagrange finite element of degree 1 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_1_1::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return true;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return false;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_1_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_1_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_1::local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_1::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_1_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_1_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_1::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = x[0][0];
    coordinates[0][1] = x[0][1];
    coordinates[0][2] = x[0][2];
    coordinates[1][0] = x[1][0];
    coordinates[1][1] = x[1][1];
    coordinates[1][2] = x[1][2];
    coordinates[2][0] = x[2][0];
    coordinates[2][1] = x[2][1];
    coordinates[2][2] = x[2][2];
    coordinates[3][0] = x[3][0];
    coordinates[3][1] = x[3][1];
    coordinates[3][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_1_1::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_1_1();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_1_2::UFC_NSEDualGradient3DLinearForm_dof_map_1_2() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_1_2::~UFC_NSEDualGradient3DLinearForm_dof_map_1_2()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_1_2::signature() const
{
    return "FFC dof map for Lagrange finite element of degree 1 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_1_2::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return true;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return false;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_1_2::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_2::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_1_2::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_2::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_2::local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_2::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_2::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_2::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_2::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_1_2::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_1_2::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1_2::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = x[0][0];
    coordinates[0][1] = x[0][1];
    coordinates[0][2] = x[0][2];
    coordinates[1][0] = x[1][0];
    coordinates[1][1] = x[1][1];
    coordinates[1][2] = x[1][2];
    coordinates[2][0] = x[2][0];
    coordinates[2][1] = x[2][1];
    coordinates[2][2] = x[2][2];
    coordinates[3][0] = x[3][0];
    coordinates[3][1] = x[3][1];
    coordinates[3][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1_2::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_1_2::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_1_2();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_1::UFC_NSEDualGradient3DLinearForm_dof_map_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_1::~UFC_NSEDualGradient3DLinearForm_dof_map_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_1::signature() const
{
    return "FFC dof map for Mixed finite element: [Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron, Lagrange finite element of degree 1 on a tetrahedron]";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_1::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return true;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return false;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1::local_dimension() const
{
    return 12;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1::num_facet_dofs() const
{
    return 9;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
    unsigned int offset = m.num_entities[0];
    dofs[4] = offset + c.entity_indices[0][0];
    dofs[5] = offset + c.entity_indices[0][1];
    dofs[6] = offset + c.entity_indices[0][2];
    dofs[7] = offset + c.entity_indices[0][3];
    offset = offset + m.num_entities[0];
    dofs[8] = offset + c.entity_indices[0][0];
    dofs[9] = offset + c.entity_indices[0][1];
    dofs[10] = offset + c.entity_indices[0][2];
    dofs[11] = offset + c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 5;
      dofs[4] = 6;
      dofs[5] = 7;
      dofs[6] = 9;
      dofs[7] = 10;
      dofs[8] = 11;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 6;
      dofs[5] = 7;
      dofs[6] = 8;
      dofs[7] = 10;
      dofs[8] = 11;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 7;
      dofs[6] = 8;
      dofs[7] = 9;
      dofs[8] = 11;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 6;
      dofs[6] = 8;
      dofs[7] = 9;
      dofs[8] = 10;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_1::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = x[0][0];
    coordinates[0][1] = x[0][1];
    coordinates[0][2] = x[0][2];
    coordinates[1][0] = x[1][0];
    coordinates[1][1] = x[1][1];
    coordinates[1][2] = x[1][2];
    coordinates[2][0] = x[2][0];
    coordinates[2][1] = x[2][1];
    coordinates[2][2] = x[2][2];
    coordinates[3][0] = x[3][0];
    coordinates[3][1] = x[3][1];
    coordinates[3][2] = x[3][2];
    coordinates[4][0] = x[0][0];
    coordinates[4][1] = x[0][1];
    coordinates[4][2] = x[0][2];
    coordinates[5][0] = x[1][0];
    coordinates[5][1] = x[1][1];
    coordinates[5][2] = x[1][2];
    coordinates[6][0] = x[2][0];
    coordinates[6][1] = x[2][1];
    coordinates[6][2] = x[2][2];
    coordinates[7][0] = x[3][0];
    coordinates[7][1] = x[3][1];
    coordinates[7][2] = x[3][2];
    coordinates[8][0] = x[0][0];
    coordinates[8][1] = x[0][1];
    coordinates[8][2] = x[0][2];
    coordinates[9][0] = x[1][0];
    coordinates[9][1] = x[1][1];
    coordinates[9][2] = x[1][2];
    coordinates[10][0] = x[2][0];
    coordinates[10][1] = x[2][1];
    coordinates[10][2] = x[2][2];
    coordinates[11][0] = x[3][0];
    coordinates[11][1] = x[3][1];
    coordinates[11][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_1::num_sub_dof_maps() const
{
    return 3;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_1::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_1_0();
      break;
    case 1:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_1_1();
      break;
    case 2:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_1_2();
      break;
    }
    return 0;
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_dof_map_2::UFC_NSEDualGradient3DLinearForm_dof_map_2() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_dof_map_2::~UFC_NSEDualGradient3DLinearForm_dof_map_2()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* UFC_NSEDualGradient3DLinearForm_dof_map_2::signature() const
{
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
}

/// Return true iff mesh entities of topological dimension d are needed
bool UFC_NSEDualGradient3DLinearForm_dof_map_2::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool UFC_NSEDualGradient3DLinearForm_dof_map_2::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void UFC_NSEDualGradient3DLinearForm_dof_map_2::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void UFC_NSEDualGradient3DLinearForm_dof_map_2::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_2::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_2::local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_2::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_2::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_2::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_2::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void UFC_NSEDualGradient3DLinearForm_dof_map_2::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    case 3:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void UFC_NSEDualGradient3DLinearForm_dof_map_2::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void UFC_NSEDualGradient3DLinearForm_dof_map_2::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int UFC_NSEDualGradient3DLinearForm_dof_map_2::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* UFC_NSEDualGradient3DLinearForm_dof_map_2::create_sub_dof_map(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_dof_map_2();
}


/// Constructor
UFC_NSEDualGradient3DLinearForm_cell_integral_0::UFC_NSEDualGradient3DLinearForm_cell_integral_0() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm_cell_integral_0::~UFC_NSEDualGradient3DLinearForm_cell_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void UFC_NSEDualGradient3DLinearForm_cell_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_02 = x[3][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    const double J_12 = x[3][1] - x[0][1];
    const double J_20 = x[1][2] - x[0][2];
    const double J_21 = x[2][2] - x[0][2];
    const double J_22 = x[3][2] - x[0][2];
      
    // Compute sub determinants
    const double d_00 = J_11*J_22 - J_12*J_21;
    const double d_01 = J_12*J_20 - J_10*J_22;
    const double d_02 = J_10*J_21 - J_11*J_20;
    
    const double d_10 = J_02*J_21 - J_01*J_22;
    const double d_11 = J_00*J_22 - J_02*J_20;
    const double d_12 = J_01*J_20 - J_00*J_21;
    
    const double d_20 = J_01*J_12 - J_02*J_11;
    const double d_21 = J_02*J_10 - J_00*J_12;
    const double d_22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d_00 + J_10*d_10 + J_20*d_20;
      
    // Compute inverse of Jacobian
    const double Jinv_00 = d_00 / detJ;
    const double Jinv_01 = d_10 / detJ;
    const double Jinv_02 = d_20 / detJ;
    const double Jinv_10 = d_01 / detJ;
    const double Jinv_11 = d_11 / detJ;
    const double Jinv_12 = d_21 / detJ;
    const double Jinv_20 = d_02 / detJ;
    const double Jinv_21 = d_12 / detJ;
    const double Jinv_22 = d_22 / detJ;
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    // Number of operations to compute element tensor = 261
    // Compute coefficients
    const double c0_0_0_0 = w[0][0];
    const double c0_0_0_1 = w[0][1];
    const double c0_0_0_2 = w[0][2];
    const double c0_0_0_3 = w[0][3];
    const double c0_0_0_4 = w[0][4];
    const double c0_0_0_5 = w[0][5];
    const double c0_0_0_6 = w[0][6];
    const double c0_0_0_7 = w[0][7];
    const double c0_0_0_8 = w[0][8];
    const double c0_0_0_9 = w[0][9];
    const double c0_0_0_10 = w[0][10];
    const double c0_0_0_11 = w[0][11];
    const double c1_0_1_0 = w[1][0];
    const double c0_1_0_0 = w[0][0];
    const double c0_1_0_1 = w[0][1];
    const double c0_1_0_2 = w[0][2];
    const double c0_1_0_3 = w[0][3];
    const double c1_1_1_0 = w[1][0];
    const double c0_2_0_0 = w[0][0];
    const double c0_2_0_1 = w[0][1];
    const double c0_2_0_2 = w[0][2];
    const double c0_2_0_3 = w[0][3];
    const double c1_2_1_0 = w[1][0];
    const double c0_3_0_4 = w[0][4];
    const double c0_3_0_5 = w[0][5];
    const double c0_3_0_6 = w[0][6];
    const double c0_3_0_7 = w[0][7];
    const double c1_3_1_0 = w[1][0];
    const double c0_4_0_4 = w[0][4];
    const double c0_4_0_5 = w[0][5];
    const double c0_4_0_6 = w[0][6];
    const double c0_4_0_7 = w[0][7];
    const double c1_4_1_0 = w[1][0];
    const double c0_5_0_8 = w[0][8];
    const double c0_5_0_9 = w[0][9];
    const double c0_5_0_10 = w[0][10];
    const double c0_5_0_11 = w[0][11];
    const double c1_5_1_0 = w[1][0];
    const double c0_6_0_8 = w[0][8];
    const double c0_6_0_9 = w[0][9];
    const double c0_6_0_10 = w[0][10];
    const double c0_6_0_11 = w[0][11];
    const double c1_6_1_0 = w[1][0];
    
    // Compute geometry tensors
    // Number of operations to compute decalrations = 162
    const double G0_0_0_0 = det*c0_0_0_0*c1_0_1_0*Jinv_00;
    const double G0_0_1_0 = det*c0_0_0_0*c1_0_1_0*Jinv_10;
    const double G0_0_2_0 = det*c0_0_0_0*c1_0_1_0*Jinv_20;
    const double G0_1_0_0 = det*c0_0_0_1*c1_0_1_0*Jinv_00;
    const double G0_2_1_0 = det*c0_0_0_2*c1_0_1_0*Jinv_10;
    const double G0_3_2_0 = det*c0_0_0_3*c1_0_1_0*Jinv_20;
    const double G0_4_0_0 = det*c0_0_0_4*c1_0_1_0*Jinv_00;
    const double G0_4_1_0 = det*c0_0_0_4*c1_0_1_0*Jinv_10;
    const double G0_4_2_0 = det*c0_0_0_4*c1_0_1_0*Jinv_20;
    const double G0_5_0_0 = det*c0_0_0_5*c1_0_1_0*Jinv_00;
    const double G0_6_1_0 = det*c0_0_0_6*c1_0_1_0*Jinv_10;
    const double G0_7_2_0 = det*c0_0_0_7*c1_0_1_0*Jinv_20;
    const double G0_8_0_0 = det*c0_0_0_8*c1_0_1_0*Jinv_00;
    const double G0_8_1_0 = det*c0_0_0_8*c1_0_1_0*Jinv_10;
    const double G0_8_2_0 = det*c0_0_0_8*c1_0_1_0*Jinv_20;
    const double G0_9_0_0 = det*c0_0_0_9*c1_0_1_0*Jinv_00;
    const double G0_10_1_0 = det*c0_0_0_10*c1_0_1_0*Jinv_10;
    const double G0_11_2_0 = det*c0_0_0_11*c1_0_1_0*Jinv_20;
    const double G1_0_0_0 = det*c0_1_0_0*c1_1_1_0*Jinv_01;
    const double G1_0_1_0 = det*c0_1_0_0*c1_1_1_0*Jinv_11;
    const double G1_0_2_0 = det*c0_1_0_0*c1_1_1_0*Jinv_21;
    const double G1_1_0_0 = det*c0_1_0_1*c1_1_1_0*Jinv_01;
    const double G1_2_1_0 = det*c0_1_0_2*c1_1_1_0*Jinv_11;
    const double G1_3_2_0 = det*c0_1_0_3*c1_1_1_0*Jinv_21;
    const double G2_0_0_0 = det*c0_2_0_0*c1_2_1_0*Jinv_02;
    const double G2_0_1_0 = det*c0_2_0_0*c1_2_1_0*Jinv_12;
    const double G2_0_2_0 = det*c0_2_0_0*c1_2_1_0*Jinv_22;
    const double G2_1_0_0 = det*c0_2_0_1*c1_2_1_0*Jinv_02;
    const double G2_2_1_0 = det*c0_2_0_2*c1_2_1_0*Jinv_12;
    const double G2_3_2_0 = det*c0_2_0_3*c1_2_1_0*Jinv_22;
    const double G3_4_0_0 = det*c0_3_0_4*c1_3_1_0*Jinv_01;
    const double G3_4_1_0 = det*c0_3_0_4*c1_3_1_0*Jinv_11;
    const double G3_4_2_0 = det*c0_3_0_4*c1_3_1_0*Jinv_21;
    const double G3_5_0_0 = det*c0_3_0_5*c1_3_1_0*Jinv_01;
    const double G3_6_1_0 = det*c0_3_0_6*c1_3_1_0*Jinv_11;
    const double G3_7_2_0 = det*c0_3_0_7*c1_3_1_0*Jinv_21;
    const double G4_4_0_0 = det*c0_4_0_4*c1_4_1_0*Jinv_02;
    const double G4_4_1_0 = det*c0_4_0_4*c1_4_1_0*Jinv_12;
    const double G4_4_2_0 = det*c0_4_0_4*c1_4_1_0*Jinv_22;
    const double G4_5_0_0 = det*c0_4_0_5*c1_4_1_0*Jinv_02;
    const double G4_6_1_0 = det*c0_4_0_6*c1_4_1_0*Jinv_12;
    const double G4_7_2_0 = det*c0_4_0_7*c1_4_1_0*Jinv_22;
    const double G5_8_0_0 = det*c0_5_0_8*c1_5_1_0*Jinv_01;
    const double G5_8_1_0 = det*c0_5_0_8*c1_5_1_0*Jinv_11;
    const double G5_8_2_0 = det*c0_5_0_8*c1_5_1_0*Jinv_21;
    const double G5_9_0_0 = det*c0_5_0_9*c1_5_1_0*Jinv_01;
    const double G5_10_1_0 = det*c0_5_0_10*c1_5_1_0*Jinv_11;
    const double G5_11_2_0 = det*c0_5_0_11*c1_5_1_0*Jinv_21;
    const double G6_8_0_0 = det*c0_6_0_8*c1_6_1_0*Jinv_02;
    const double G6_8_1_0 = det*c0_6_0_8*c1_6_1_0*Jinv_12;
    const double G6_8_2_0 = det*c0_6_0_8*c1_6_1_0*Jinv_22;
    const double G6_9_0_0 = det*c0_6_0_9*c1_6_1_0*Jinv_02;
    const double G6_10_1_0 = det*c0_6_0_10*c1_6_1_0*Jinv_12;
    const double G6_11_2_0 = det*c0_6_0_11*c1_6_1_0*Jinv_22;
    
    // Compute element tensor
    // Number of operations to compute tensor = 99
    A[0] = -0.166666666666666*G0_0_0_0 - 0.166666666666666*G0_0_1_0 - 0.166666666666666*G0_0_2_0 + 0.166666666666666*G0_1_0_0 + 0.166666666666666*G0_2_1_0 + 0.166666666666666*G0_3_2_0;
    A[1] = -0.166666666666666*G0_4_0_0 - 0.166666666666666*G0_4_1_0 - 0.166666666666666*G0_4_2_0 + 0.166666666666666*G0_5_0_0 + 0.166666666666666*G0_6_1_0 + 0.166666666666666*G0_7_2_0;
    A[2] = -0.166666666666666*G0_8_0_0 - 0.166666666666666*G0_8_1_0 - 0.166666666666666*G0_8_2_0 + 0.166666666666666*G0_9_0_0 + 0.166666666666666*G0_10_1_0 + 0.166666666666666*G0_11_2_0;
    A[3] = -0.166666666666666*G1_0_0_0 - 0.166666666666666*G1_0_1_0 - 0.166666666666666*G1_0_2_0 + 0.166666666666666*G1_1_0_0 + 0.166666666666666*G1_2_1_0 + 0.166666666666666*G1_3_2_0;
    A[4] = -0.166666666666666*G3_4_0_0 - 0.166666666666666*G3_4_1_0 - 0.166666666666666*G3_4_2_0 + 0.166666666666666*G3_5_0_0 + 0.166666666666666*G3_6_1_0 + 0.166666666666666*G3_7_2_0;
    A[5] = -0.166666666666666*G5_8_0_0 - 0.166666666666666*G5_8_1_0 - 0.166666666666666*G5_8_2_0 + 0.166666666666666*G5_9_0_0 + 0.166666666666666*G5_10_1_0 + 0.166666666666666*G5_11_2_0;
    A[6] = -0.166666666666666*G2_0_0_0 - 0.166666666666666*G2_0_1_0 - 0.166666666666666*G2_0_2_0 + 0.166666666666666*G2_1_0_0 + 0.166666666666666*G2_2_1_0 + 0.166666666666666*G2_3_2_0;
    A[7] = -0.166666666666666*G4_4_0_0 - 0.166666666666666*G4_4_1_0 - 0.166666666666666*G4_4_2_0 + 0.166666666666666*G4_5_0_0 + 0.166666666666666*G4_6_1_0 + 0.166666666666666*G4_7_2_0;
    A[8] = -0.166666666666666*G6_8_0_0 - 0.166666666666666*G6_8_1_0 - 0.166666666666666*G6_8_2_0 + 0.166666666666666*G6_9_0_0 + 0.166666666666666*G6_10_1_0 + 0.166666666666666*G6_11_2_0;
}

/// Constructor
UFC_NSEDualGradient3DLinearForm::UFC_NSEDualGradient3DLinearForm() : ufc::form()
{
    // Do nothing
}

/// Destructor
UFC_NSEDualGradient3DLinearForm::~UFC_NSEDualGradient3DLinearForm()
{
    // Do nothing
}

/// Return a string identifying the form
const char* UFC_NSEDualGradient3DLinearForm::signature() const
{
    return "w0_a0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]w1_a2[0](dXa1[0, 1, 2]/dx0) | ((d/dXa1[0, 1, 2])va0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11][b0[0, 1, 2]])*vi0[0, 1, 2, 3, 4, 5, 6, 7, 8][b0[0, 1, 2]]*va2[0]*dX(0) + w0_a0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]w1_a2[0](dXa1[0, 1, 2]/dx1) | ((d/dXa1[0, 1, 2])va0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11][0])*vi0[0, 1, 2, 3, 4, 5, 6, 7, 8][3]*va2[0]*dX(0) + w0_a0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]w1_a2[0](dXa1[0, 1, 2]/dx2) | ((d/dXa1[0, 1, 2])va0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11][0])*vi0[0, 1, 2, 3, 4, 5, 6, 7, 8][6]*va2[0]*dX(0) + w0_a0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]w1_a2[0](dXa1[0, 1, 2]/dx1) | ((d/dXa1[0, 1, 2])va0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11][1])*vi0[0, 1, 2, 3, 4, 5, 6, 7, 8][4]*va2[0]*dX(0) + w0_a0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]w1_a2[0](dXa1[0, 1, 2]/dx2) | ((d/dXa1[0, 1, 2])va0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11][1])*vi0[0, 1, 2, 3, 4, 5, 6, 7, 8][7]*va2[0]*dX(0) + w0_a0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]w1_a2[0](dXa1[0, 1, 2]/dx1) | ((d/dXa1[0, 1, 2])va0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11][2])*vi0[0, 1, 2, 3, 4, 5, 6, 7, 8][5]*va2[0]*dX(0) + w0_a0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]w1_a2[0](dXa1[0, 1, 2]/dx2) | ((d/dXa1[0, 1, 2])va0[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11][2])*vi0[0, 1, 2, 3, 4, 5, 6, 7, 8][8]*va2[0]*dX(0)";
}

/// Return the rank of the global tensor (r)
unsigned int UFC_NSEDualGradient3DLinearForm::rank() const
{
    return 1;
}

/// Return the number of coefficients (n)
unsigned int UFC_NSEDualGradient3DLinearForm::num_coefficients() const
{
    return 2;
}

/// Return the number of cell integrals
unsigned int UFC_NSEDualGradient3DLinearForm::num_cell_integrals() const
{
    return 1;
}
  
/// Return the number of exterior facet integrals
unsigned int UFC_NSEDualGradient3DLinearForm::num_exterior_facet_integrals() const
{
    return 0;
}
  
/// Return the number of interior facet integrals
unsigned int UFC_NSEDualGradient3DLinearForm::num_interior_facet_integrals() const
{
    return 0;
}

/// Create a new finite element for argument function i
ufc::finite_element* UFC_NSEDualGradient3DLinearForm::create_finite_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_0();
      break;
    case 1:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_1();
      break;
    case 2:
      return new UFC_NSEDualGradient3DLinearForm_finite_element_2();
      break;
    }
    return 0;
}

/// Create a new dof map for argument function i
ufc::dof_map* UFC_NSEDualGradient3DLinearForm::create_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_0();
      break;
    case 1:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_1();
      break;
    case 2:
      return new UFC_NSEDualGradient3DLinearForm_dof_map_2();
      break;
    }
    return 0;
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* UFC_NSEDualGradient3DLinearForm::create_cell_integral(unsigned int i) const
{
    return new UFC_NSEDualGradient3DLinearForm_cell_integral_0();
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* UFC_NSEDualGradient3DLinearForm::create_exterior_facet_integral(unsigned int i) const
{
    return 0;
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* UFC_NSEDualGradient3DLinearForm::create_interior_facet_integral(unsigned int i) const
{
    return 0;
}

