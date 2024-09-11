/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this;
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <iostream>
#include "../../indexSquaredDistance.h"
#include "../../../Common/Containers/stack.h"

namespace Ponca {
template <typename Traits> class KdTreeBase;

template <typename Traits>
class KdTreeQuery
{
public:
    using DataPoint  = typename Traits::DataPoint;
    using IndexType  = typename Traits::IndexType;
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;

    explicit inline KdTreeQuery(const KdTreeBase<Traits>* kdtree) : m_kdtree( kdtree ), m_stack() {}

protected:
    /// \brief Init stack for a new search
    inline void reset() {
        m_stack.clear();
        m_stack.push({0,0});
    }

    /// [KdTreeQuery kdtree type]
    const KdTreeBase<Traits>* m_kdtree { nullptr };
    /// [KdTreeQuery kdtree type]
    Stack<IndexSquaredDistance<IndexType, Scalar>, 2 * Traits::MAX_DEPTH> m_stack;

    /// \return false if the kdtree is empty
    template<typename LeafPreparationFunctor,
            typename DescentDistanceThresholdFunctor,
            typename SkipIndexFunctor,
            typename ProcessNeighborFunctor>
    bool search_internal(const VectorType& point,
                         LeafPreparationFunctor prepareLeafTraversal,
                         DescentDistanceThresholdFunctor descentDistanceThreshold,
                         SkipIndexFunctor skipFunctor,
                         ProcessNeighborFunctor processNeighborFunctor
                         )
    {
        const auto& nodes  = m_kdtree->nodes();
        const auto& points = m_kdtree->points();
        // if(DataPoint::Dim==5){
        //   std::cout<<"search_internal "<<point<<std::endl;
        //   std::cout<<"stack size "<<m_stack.size()<<std::endl;
        //   std::cout<<"nodes size "<<nodes.size()<<std::endl;
        //   std::cout<<"points size "<<points.size()<<std::endl;
        // }

        if (nodes.empty() || points.empty() || m_kdtree->sample_count() == 0){
        // if(DataPoint::Dim==5){
        //   std::cout<<"emptuy, exit"<<std::endl;
        // }
            return false;
        }

        while(!m_stack.empty())
        {
            auto& qnode = m_stack.top();
            // if(DataPoint::Dim==5){
            //   std::cout<<" qnode "<<qnode.index<<" "<<qnode.squared_distance<<std::endl;
            // }
            const auto& node = nodes[qnode.index];
            // if(DataPoint::Dim==5){
            // std::cout<<" corres node "<<qnode.index<<" "<<node.is_leaf()<<" "<<node.inner_first_child_id()<<" "<<node.inner_split_dim()<<" "<<node.inner_split_value()<<std::endl;
            // }

            if(qnode.squared_distance < descentDistanceThreshold())
            {
                if(node.is_leaf())
                {
                    // if(DataPoint::Dim==5){
                    // std::cout<<" corres node is leaf"<<std::endl;
                    // }
                    m_stack.pop();
                    IndexType start = node.leaf_start();
                    IndexType end = node.leaf_start() + node.leaf_size();
                    prepareLeafTraversal(start, end);
                    for(IndexType i=start; i<end; ++i)
                    {
                      // if(DataPoint::Dim==5){
                      //   std::cout<<" sample "<<i<<std::endl;
                      // }
                        IndexType idx = m_kdtree->pointFromSample(i);
        // if(DataPoint::Dim==5){
        //   std::cout<<"point from sample "<<i<<" -> "<<idx<<std::endl;
        // }
                        if(skipFunctor(idx)) continue;

                        Scalar d = (point - points[idx].pos()).squaredNorm();

                        if(d < descentDistanceThreshold())
                        {
                            if( processNeighborFunctor( idx, i, d )) return false;
                        }
                    }
                }
                else
                {
                    // replace the stack top by the farthest and push the closest
                    Scalar newOff = point[node.inner_split_dim()] - node.inner_split_value();
                    // if(DataPoint::Dim==5){
                    //   std::cout<<"newOff "<<node.inner_split_dim()<<" "<<node.inner_split_value()<<" "<<newOff<<std::endl;
                    // }
                    m_stack.push();
                    if(newOff < 0)
                    {
                        m_stack.top().index = node.inner_first_child_id();
                        qnode.index         = node.inner_first_child_id()+1;
                    }
                    else
                    {
                        m_stack.top().index = node.inner_first_child_id()+1;
                        qnode.index         = node.inner_first_child_id();
                    }
                    m_stack.top().squared_distance = qnode.squared_distance;
                    qnode.squared_distance         = newOff*newOff;
                    // if(DataPoint::Dim==5){
                    //   std::cout<<"m_stack push "<<m_stack.top().index<<" "<<m_stack.top().squared_distance<<std::endl;
                    //   std::cout<<"qnode is now "<<qnode.index<<" "<<qnode.squared_distance<<std::endl;
                    // }
                }
            }
            else
            {
                m_stack.pop();
            }
        }
        // if(DataPoint::Dim==5){
        //   std::cout<<"search_internal done"<<std::endl;
        // }
        return true;
    }
};
} // namespace Ponca
