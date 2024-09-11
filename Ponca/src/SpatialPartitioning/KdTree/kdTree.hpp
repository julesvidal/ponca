/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// KdTree ----------------------------------------------------------------------

template<typename Traits>
template<typename PointUserContainer, typename Converter>
inline void KdTreeBase<Traits>::build(PointUserContainer&& points, Converter c)
{
    IndexContainer ids(points.size());
    std::iota(ids.begin(), ids.end(), 0);
    this->buildWithSampling(std::forward<PointUserContainer>(points), std::move(ids), std::move(c));
}

template<typename Traits>
void KdTreeBase<Traits>::clear()
{
    m_points.clear();
    m_nodes.clear();
    m_indices.clear();
    m_leaf_count = 0;
}

template<typename Traits>
bool KdTreeBase<Traits>::valid() const
{
    if (m_points.empty())
        return m_nodes.empty() && m_indices.empty();

    if(m_nodes.empty() || m_indices.empty())
    {
        return false;
    }

    std::vector<bool> b(point_count(), false);
    for(IndexType idx : m_indices)
    {
        if(idx < 0 || point_count() <= idx || b[idx])
        {
            return false;
        }
        b[idx] = true;
    }

    for(NodeIndexType n=0;n<node_count();++n)
    {
        const NodeType& node = m_nodes[n];
        if(node.is_leaf())
        {
            if(sample_count() <= node.leaf_start() || node.leaf_start()+node.leaf_size() > sample_count())
            {
                return false;
            }
        }
        else
        {
            if(node.inner_split_dim() < 0 || DataPoint::Dim-1 < node.inner_split_dim())
            {
                return false;
            }
            if(node_count() <= node.inner_first_child_id() || node_count() <= node.inner_first_child_id()+1)
            {
                return false;
            }
        }
    }

    return true;
}

template<typename Traits>
void KdTreeBase<Traits>::print(std::ostream& os, bool verbose) const
{
    os << "KdTree:";
    os << "\n  MaxNodes: " << MAX_NODE_COUNT;
    os << "\n  MaxPoints: " << MAX_POINT_COUNT;
    os << "\n  MaxDepth: " << MAX_DEPTH;
    os << "\n  PointCount: " << point_count();
    os << "\n  SampleCount: " << sample_count();
    os << "\n  NodeCount: " << node_count();

    if (!verbose)
    {
        return;
    }

    os << "\n  Samples: [";
    static constexpr IndexType SAMPLES_PER_LINE = 10;
    for (IndexType i = 0; i < sample_count(); ++i)
    {
        os << (i == 0 ? "" : ",");
        os << (i % SAMPLES_PER_LINE == 0 ? "\n    " : " ");
        os << m_indices[i];
    }

    os << "]\n  Nodes:";
    for (NodeIndexType n = 0; n < node_count(); ++n)
    {
        const NodeType& node = m_nodes[n];
        if (node.is_leaf())
        {
            os << "\n    - Type: Leaf";
            os << "\n      Start: " << node.leaf_start();
            os << "\n      Size: " << node.leaf_size();
        }
        else
        {
            os << "\n    - Type: Inner";
            os << "\n      SplitDim: " << node.inner_split_dim();
            os << "\n      SplitValue: " << node.inner_split_value();
            os << "\n      FirstChild: " << node.inner_first_child_id();
        }
    }
}

template<typename Traits>
template<typename PointUserContainer, typename IndexUserContainer, typename Converter>
inline void KdTreeBase<Traits>::buildWithSampling(PointUserContainer&& points,
                                                  IndexUserContainer sampling,
                                                  Converter c)
{
    PONCA_DEBUG_ASSERT(points.size() <= MAX_POINT_COUNT);
    this->clear();

    // Move, copy or convert input samples
    c(std::forward<PointUserContainer>(points), m_points);

    m_nodes = NodeContainer();
    m_nodes.reserve(4 * point_count() / m_min_cell_size);

        // if(DataPoint::Dim==5){
        //   std::cout<<"m_nodes reserve "<<4 * point_count() / m_min_cell_size<<" "<<point_count()<<" "<<m_min_cell_size<<std::endl;
        // }
    m_nodes.emplace_back();

    m_indices = std::move(sampling);

    
    // if(DataPoint::Dim==5){
    //   std::cout<<"build with sampling"<<std::endl;
    //   std::cout<<"sample count"<<sample_count()<<std::endl;
    //   std::cout<<"point count"<<point_count()<<std::endl;
    // }

    this->build_rec(0, 0, sample_count(), 1);

    if(DataPoint::Dim==5){
      // std::cout<<"kd tree check nodes"<<std::endl;
      // std::cout<<m_nodes.size()<<" nodes"<<std::endl;
      for(int inode=0; inode<m_nodes.size(); inode++){
        auto& node = m_nodes[inode];
        // std::cout<<" node "<<inode<<" "<<node.is_leaf()<<" "<<node.inner_first_child_id() \
          <<" "<<node.inner_split_value()<<" "<<node.inner_split_dim()<<std::endl;
      }
    }
    PONCA_DEBUG_ASSERT(this->valid());

}

template<typename Traits>
void KdTreeBase<Traits>::build_rec(NodeIndexType node_id, IndexType start, IndexType end, int level)
{
            if(DataPoint::Dim==5){
  // std::cout<<" node id "<<node_id<<" "<<&m_nodes[node_id]<<std::endl;
            }
  if(node_id>m_nodes.size() or node_id<0){
    std::cout<<"ERROR NODE ID "<<node_id<<" outta "<<m_nodes.size()<<std::endl;
  }
    NodeType& node = m_nodes[node_id];
    AabbType aabb;
    for(IndexType i=start; i<end; ++i){
        aabb.extend(m_points[m_indices[i]].pos());
            if(DataPoint::Dim==5){
              // std::cout<<" aabb "<<aabb.dim()<<" "<<aabb.volume()<<std::endl;
            }
    }

    m_nodes[node_id].set_is_leaf(
        end-start <= m_min_cell_size ||
        level >= Traits::MAX_DEPTH ||
        // Since we add 2 nodes per inner node we need to stop if we can't add
        // them both
        (NodeIndexType)m_nodes.size() > MAX_NODE_COUNT - 2 ||
        aabb.min().isApprox(aabb.max()));

    m_nodes[node_id].configure_range(start, end-start, aabb);
            // if(DataPoint::Dim==5){
          // std::cout<<"  build node "<<std::endl;
        // }

    if (m_nodes[node_id].is_leaf())
    {
            // if(DataPoint::Dim==5){
            //   auto& node2 = m_nodes[node_id];
          // std::cout<<"  build leaf node "<<node_id<<" "<<node2.is_leaf()<<" "<<node2.inner_first_child_id()<<" "<<node2.inner_split_value()<<" "<<node2.inner_split_dim()<<std::endl;
        // }
        ++m_leaf_count;
    }
    else
    {
        int split_dim = 0;
        (Scalar(0.5) * aabb.diagonal()).maxCoeff(&split_dim);
        m_nodes[node_id].configure_inner(aabb.center()[split_dim], m_nodes.size(), split_dim);

        auto first_child_id = m_nodes[node_id].inner_first_child_id();
        auto split_value = m_nodes[node_id].inner_split_value();
            // if(DataPoint::Dim==5){
            //   auto& node2 = m_nodes[node_id];
          // std::cout<<"  build inner node "<<node_id<<" "<<node2.is_leaf()<<" "<<node2.inner_first_child_id()<<" "<<node2.inner_split_value()<<" "<<node2.inner_split_dim()<<std::endl;
          // std::cout<<"  build inner node2 "<<node_id<<" "<<node.is_leaf()<<" "<<node.inner_first_child_id()<<" "<<node.inner_split_value()<<" "<<node.inner_split_dim()<<std::endl;
          // std::cout<<"  build inner node3 "<<node_id<<" "<<m_nodes[node_id].is_leaf()<<" "<<m_nodes[node_id].inner_first_child_id()<<" "<<m_nodes[node_id].inner_split_value()<<" "<<m_nodes[node_id].inner_split_dim()<<std::endl;
          // std::cout<<" node id "<<node_id<<" adress "<<&m_nodes[node_id]<<std::endl;
          // std::cout<<" aabb "<<aabb.diagonal().norm()<<std::endl;
        // }
        // if(DataPoint::Dim==5){
        // std::cout<<" size-cap "<<m_nodes.capacity()-m_nodes.size()<<std::endl;
        // std::cout<<"m_nodes size"<<node_count()<<std::endl;
        // }
        m_nodes.emplace_back();
        m_nodes.emplace_back();

            // if(DataPoint::Dim==5){
            //   auto& node2 = m_nodes[node_id];
          // std::cout<<"       inner node "<<node_id<<" has changed to  "<<node2.is_leaf()<<" "<<node2.inner_first_child_id()<<" "<<node2.inner_split_value()<<" "<<node2.inner_split_dim()<<std::endl;
          // std::cout<<"       inner node2 "<<node_id<<" has changed to  "<<node.is_leaf()<<" "<<node.inner_first_child_id()<<" "<<node.inner_split_value()<<" "<<node.inner_split_dim()<<std::endl;
          // std::cout<<"       inner node3 "<<node_id<<" has changed to  "<<m_nodes[node_id].is_leaf()<<" "<<m_nodes[node_id].inner_first_child_id()<<" "<<m_nodes[node_id].inner_split_value()<<" "<<m_nodes[node_id].inner_split_dim()<<std::endl;
          // std::cout<<" node id "<<node_id<<" adress "<<&m_nodes[node_id]<<std::endl;
        // }

        IndexType mid_id = this->partition(start, end, split_dim, m_nodes[node_id].inner_split_value());
        // if(DataPoint::Dim==5){
        //   std::cout<<"  mid_id "<<mid_id<<std::endl;
        //   auto& node2 = m_nodes[node_id];
        //   std::cout<<"       inner node "<<node_id<<" has now changed to  "<<node2.is_leaf()<<" "<<node2.inner_first_child_id()<<" "<<node2.inner_split_value()<<" "<<node2.inner_split_dim()<<std::endl;
        //   std::cout<<"       inner node2 "<<node_id<<" has now changed to  "<<node.is_leaf()<<" "<<node.inner_first_child_id()<<" "<<node.inner_split_value()<<" "<<node.inner_split_dim()<<std::endl;
        //   std::cout<<"       inner node3 "<<node_id<<" has now changed to  "<<m_nodes[node_id].is_leaf()<<" "<<m_nodes[node_id].inner_first_child_id()<<" "<<m_nodes[node_id].inner_split_value()<<" "<<m_nodes[node_id].inner_split_dim()<<std::endl;
        //   std::cout<<" node id "<<node_id<<" adress "<<&m_nodes[node_id]<<std::endl;
        // }
        // if(DataPoint::Dim==5){
        //   std::cout<<" recursive call "<<node_id<<" with "<<first_child_id<<std::endl;
        // }
        build_rec(first_child_id, start, mid_id, level+1);
        // if(DataPoint::Dim==5){
        //   std::cout<<"  mid_id "<<mid_id<<std::endl;
        //   auto& node2 = m_nodes[node_id];
        //   std::cout<<"       inner node "<<node_id<<" has now changed to  "<<node2.is_leaf()<<" "<<node2.inner_first_child_id()<<" "<<node2.inner_split_value()<<" "<<node2.inner_split_dim()<<std::endl;
        //   std::cout<<"       inner node2 "<<node_id<<" has now changed to  "<<node.is_leaf()<<" "<<node.inner_first_child_id()<<" "<<node.inner_split_value()<<" "<<node.inner_split_dim()<<std::endl;
        //   std::cout<<"       inner node3 "<<node_id<<" has now changed to  "<<m_nodes[node_id].is_leaf()<<" "<<m_nodes[node_id].inner_first_child_id()<<" "<<m_nodes[node_id].inner_split_value()<<" "<<m_nodes[node_id].inner_split_dim()<<std::endl;
        //   std::cout<<" node id "<<node_id<<" adress "<<&m_nodes[node_id]<<std::endl;
        // }
        // if(DataPoint::Dim==5){
        //   std::cout<<" recursive call "<<node_id<<" with "<<first_child_id+1<<std::endl;
        // }
        build_rec(first_child_id+1, mid_id, end, level+1);
            // if(DataPoint::Dim==5){
            //   auto& node2 = m_nodes[node_id];
          // std::cout<<"       inner node "<<node_id<<" is now  "<<node2.is_leaf()<<" "<<node2.inner_first_child_id()<<" "<<node2.inner_split_value()<<" "<<node2.inner_split_dim()<<std::endl;
          // std::cout<<"       inner node2 "<<node_id<<" is now  "<<node.is_leaf()<<" "<<node.inner_first_child_id()<<" "<<node.inner_split_value()<<" "<<node.inner_split_dim()<<std::endl;
          // std::cout<<"       inner node3 "<<node_id<<" is now  "<<m_nodes[node_id].is_leaf()<<" "<<m_nodes[node_id].inner_first_child_id()<<" "<<m_nodes[node_id].inner_split_value()<<" "<<m_nodes[node_id].inner_split_dim()<<std::endl;
          // std::cout<<" node id "<<node_id<<" adress "<<&m_nodes[node_id]<<std::endl;
        // }
        m_nodes[node_id].configure_inner(split_value, first_child_id, split_dim);
        // m_nodes[node_id].inner_first_child_id()=first_child_id;
        // m_nodes[node_id].inner_split_value()=first_child_id;
        // m_nodes[node_id].inner_split_dim()=split_dim;
            // if(DataPoint::Dim==5){
            //   auto& node2 = m_nodes[node_id];
          // std::cout<<"       inner node "<<node_id<<" is finally  "<<node2.is_leaf()<<" "<<node2.inner_first_child_id()<<" "<<node2.inner_split_value()<<" "<<node2.inner_split_dim()<<std::endl;
          // std::cout<<"       inner node2 "<<node_id<<" is finally  "<<node.is_leaf()<<" "<<node.inner_first_child_id()<<" "<<node.inner_split_value()<<" "<<node.inner_split_dim()<<std::endl;
          // std::cout<<"       inner node3 "<<node_id<<" is finally  "<<m_nodes[node_id].is_leaf()<<" "<<m_nodes[node_id].inner_first_child_id()<<" "<<m_nodes[node_id].inner_split_value()<<" "<<m_nodes[node_id].inner_split_dim()<<std::endl;
          // std::cout<<" node id "<<node_id<<" adress "<<&m_nodes[node_id]<<std::endl;
            // }
    }
    // if(DataPoint::Dim==5){
    //   std::cout<<"m_nodes size"<<node_count()<<std::endl;
    // }
}

template<typename Traits>
auto KdTreeBase<Traits>::partition(IndexType start, IndexType end, int dim, Scalar value)
    -> IndexType
{
    const auto& points = m_points;
    auto& indices  = m_indices;
    
    auto it = std::partition(indices.begin()+start, indices.begin()+end, [&](IndexType i)
    {
        return points[i].pos()[dim] < value;
    });

    auto distance = std::distance(m_indices.begin(), it);
    
    return static_cast<IndexType>(distance);
}
