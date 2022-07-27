/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "./defines.h"

#include <Eigen/Dense>

namespace Ponca
{

/*!
 * \brief Extension to compute the best fit quadric on 3d points expressed as \f$f(u,v)=h\f$
 *
 * \note This procedure requires at least two passes, the first one for plane fitting,
 * the second one for quadric fitting.
 * \warning This class is valid only in 3D.
 *
 * \note This class mixes the primitive (MongePatch) and its fitting procedure.
 *       Could makes sense to split the two
 * \ingroup fitting
 */
template < class DataPoint, class _WFunctor, typename T>
class MongePatch : public T
{
private:
    using Base = T;

protected:
    enum
    {
        Check = Base::PROVIDES_PLANE && Base::PROVIDES_TANGENT_PLANE_BASIS
    };

public:
    using Scalar     = typename Base::Scalar;     /*!< \brief Inherited scalar type*/
    using VectorType = typename Base::VectorType; /*!< \brief Inherited vector type*/
    using WFunctor   = typename Base::WFunctor;   /*!< \brief Weight Function*/

    using SampleMatrix = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    using Vector6      = Eigen::Matrix<Scalar,6,1>;

protected:
    SampleMatrix m_A; /*!< \brief Quadric input samples */
    Vector6      m_x {Vector6::Zero()};      /*!< \brief Quadric parameters */
    Vector6      m_b {Vector6::Zero()};      /*!< \brief Observations */

    bool m_planeIsReady {false};
public:
    /*! \brief Default constructor */
    PONCA_MULTIARCH inline MongePatch() = default;

    PONCA_EXPLICIT_CAST_OPERATORS(MongePatch,mongePatch)

    /**************************************************************************/
    /* Initialization                                                         */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::init() */
    PONCA_MULTIARCH inline void init (const VectorType& _evalPos);

    /**************************************************************************/
    /* Processing                                                             */
    /**************************************************************************/
    /*! \copydoc Concept::FittingProcedureConcept::addLocalNeighbor() */
    PONCA_MULTIARCH inline bool addLocalNeighbor(Scalar w, const VectorType &localQ, const DataPoint &attributes);

    /*! \copydoc Concept::FittingProcedureConcept::finalize() */
    PONCA_MULTIARCH inline FIT_RESULT finalize();


    //! \brief Returns an estimate of the mean curvature
    PONCA_MULTIARCH inline Scalar kMean() const;

    //! \brief Returns an estimate of the Gaussian curvature
    PONCA_MULTIARCH inline Scalar GaussianCurvature() const;

    PONCA_MULTIARCH inline Scalar evalUV(Scalar u, Scalar v) const {
      return h_uu()*u*u + h_vv()*v*v + h_uv()*u*v + h_u()*u + h_v()*v + h_c();
    }

    /*! \brief Value of the scalar field at the evaluation point */
    PONCA_MULTIARCH inline Scalar potential(const VectorType& _q) const {
      VectorType x = Base::worldToTangentPlane(_q);
      return evalUV(*(x.data()+1),*(x.data()+2)) - *(x.data());
    }

    //! \brief Orthogonal projecting on the patch, such that h = f(u,v)
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        VectorType x = Base::worldToTangentPlane(_q);
        *(x.data()) = evalUV(*(x.data()+1),*(x.data()+2));
        return Base::tangentPlaneToWorld(x);
    }

    PONCA_MULTIARCH inline const Scalar & h_uu () const { return *(m_x.data()); }
    PONCA_MULTIARCH inline const Scalar & h_vv () const { return *(m_x.data()+1); }
    PONCA_MULTIARCH inline const Scalar & h_uv () const { return *(m_x.data()+2); }
    PONCA_MULTIARCH inline const Scalar & h_u  () const { return *(m_x.data()+3); }
    PONCA_MULTIARCH inline const Scalar & h_v  () const { return *(m_x.data()+4); }
    PONCA_MULTIARCH inline const Scalar & h_c  () const { return *(m_x.data()+5); }

};

#include "mongePatch.hpp"

} //namespace Ponca
