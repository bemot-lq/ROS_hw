#include <graph_searcher.h>

using namespace std;
using namespace Eigen;

void gridPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);//最低
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);//最高
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;//_max_x_id = (int)(_x_size * _inv_resolution);
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;//点云数据采用一维数组，按(0,0,0),(0,0,1)存储(x,y,z),栅格地图(x,y,z)对应点云坐标在数组中的位置
    //x*GLYZ_SIZE+y*GLZ_SIZE+z

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];//点云数组
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));//给data全赋0
    
    //将不同坐标系下的同一点放在类成员中进行管理
    /*
        GridNode(Eigen::Vector3i _index, Eigen::Vector3d _coord){  
		id = 0;
		is_path = false;
		index = _index;
		coord = _coord;
		dir   = Eigen::Vector3i::Zero();

		gScore = inf;
		fScore = inf;
		cameFrom = NULL;
    
    */
    GridNodeMap = new GridNodePtr ** [GLX_SIZE];//使用指针构建三维数组，GridNodeMap只是一个三维指针数组和GridNode类型不同
    for(int i = 0; i < GLX_SIZE; i++){
        GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
        for(int j = 0; j < GLY_SIZE; j++){
            GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
            for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);//将栅格坐标对应到点云坐标
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);//前面的是栅格，后面的是点云
            }
        }
    }
}

void gridPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;//data点云一维数组，对相应位置赋1
}




//栅格坐标系均大于0,转换点云需要加上最小值，移动至中心点，乘上分辨率
inline Vector3d gridPathFinder::gridIndex2coord(const Vector3i & index) const
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}
//点云不一定在正中间，除以分辨率，与取整有关。GLX_SIZE为size除分辨率有关，点云分辨率高于栅格的
inline Vector3i gridPathFinder::coord2gridIndex(const Vector3d & pt) const
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}
//点云分辨率高于栅格
Eigen::Vector3d gridPathFinder::coordRounding(const Eigen::Vector3d & coord) const
{
    return gridIndex2coord(coord2gridIndex(coord));
}

void gridPathFinder::astarGraphSearch(const Eigen::Vector3d start_pt, const Eigen::Vector3d target_pt)
{
    Vector3d start_pt_n=coordRounding(start_pt);
    Vector3d target_pt_n=coordRounding(target_pt);
    Vector3i start_idx=coord2gridIndex(start_pt);
    Vector3i target_idx=coord2gridIndex(target_pt);

    

    //将初始点放入开集（手动设置）
    start=GridNodeMap[start_idx(0)][start_idx(1)][start_idx(2)];
    target=GridNodeMap[target_idx(0)][target_idx(1)][target_idx(2)];

    start->id=1;
    start->fScore=fabs(start->coord(0)-target->coord(0))+
                  fabs(start->coord(1)-target->coord(1))+
                  fabs(start->coord(2)-target->coord(2));
    openSet.insert({start->gScore+start->fScore,start});

    while (!openSet.empty()||notReachTarget)
    {
        //从开集中读取点，进行分支，放入开集（非闭放入开集）
        GridNodePtr node=openSet.begin()->second;
        //pop第一个
        openSet.erase(openSet.begin());
        expandNodes(node);

    }

    if (notReachTarget==false)
    {
        cout<<"Reach Target"<<endl;
    }
    else
    {
        cout<<"Reach False"<<endl;
    }
    //1.开集扩展函数2.开集读取函数
    //检测是否到达终点，没有重复上一步骤

    //到达终点，累进父节点读出路径
    showGridPathFun();
}

bool gridPathFinder::judeReachTarget(const Eigen::Vector3d pt_a,const Eigen::Vector3d pt_b)
{
    if (pt_a(0)==pt_b(0)&&pt_a(1)==pt_b(1)&&pt_a(2)==pt_b(2))
    {
        return true;
    }
    else
    {
        return false;
    }
}


void gridPathFinder::expandNodes(GridNodePtr node)
{
    Vector3i idx=node->index;
    if (idx(0)+1<GLX_SIZE)
    {
        GridNodePtr node_1=GridNodeMap[idx(0)+1][idx(1)][idx(2)];
        if (judeReachTarget(target->coord,node_1->coord))
        {
            notReachTarget=false;
            target->cameFrom=node;
            return;
        }
        if (data[(idx(0)+1)*GLYZ_SIZE+idx(1)*GLZ_SIZE+idx(2)]==0)
        {
            
            expandNode(node,node_1);
        }
        
    }
    if (idx(0)-1>=0)
    {
        GridNodePtr node_2=GridNodeMap[idx(0)-1][idx(1)][idx(2)];
        if (judeReachTarget(target->coord,node_2->coord))
        {
            notReachTarget=false;
            target->cameFrom=node;
            return;
        }
        if (data[(idx(0)-1)*GLYZ_SIZE+idx(1)*GLZ_SIZE+idx(2)]==0)
        {
            GridNodePtr node_2=GridNodeMap[idx(0)-1][idx(1)][idx(2)];
            expandNode(node,node_2);
        }
        
    }
    if (idx(1)+1<GLY_SIZE)
    {
        GridNodePtr node_3=GridNodeMap[idx(0)][idx(1)+1][idx(2)];
        if (judeReachTarget(target->coord,node_3->coord))
        {
            notReachTarget=false;
            target->cameFrom=node;
            return;
        }
        if (data[idx(0)*GLYZ_SIZE+(idx(1)+1)*GLZ_SIZE+idx(2)]==0)
        {
            expandNode(node,node_3);
        }
        
    }
    if (idx(1)-1>=0)
    {
        GridNodePtr node_4=GridNodeMap[idx(0)][idx(1)-1][idx(2)];
        if (judeReachTarget(target->coord,node_4->coord))
        {
            notReachTarget=false;
            target->cameFrom=node;
            return;
        }
        if (data[idx(0)*GLYZ_SIZE+(idx(1)-1)*GLZ_SIZE+idx(2)]==0)
        {
            expandNode(node,node_4);
        }
        
    }    
    if (idx(2)+1<GLZ_SIZE)
    {
        GridNodePtr node_5=GridNodeMap[idx(0)][idx(1)][idx(2)+1];
        if (judeReachTarget(target->coord,node_5->coord))
        {
            notReachTarget=false;
            target->cameFrom=node;
            return;
        }
        if (data[idx(0)*GLYZ_SIZE+idx(1)*GLZ_SIZE+(idx(2)+1)]==0)
        {
            expandNode(node,node_5);
        }
        
    }
    if (idx(2)-1>=0)
    {
        GridNodePtr node_6=GridNodeMap[idx(0)][idx(1)][idx(2)-1];
        if (judeReachTarget(target->coord,node_6->coord))
        {
            notReachTarget=false;
            target->cameFrom=node;
            return;
        }
        if (data[idx(0)*GLYZ_SIZE+idx(1)*GLZ_SIZE+(idx(2)-1)]==0)
        {
            GridNodePtr node_6=GridNodeMap[idx(0)][idx(1)][idx(2)-1];
            expandNode(node,node_6);
        }
        
    }
    

}

void gridPathFinder::expandNode(GridNodePtr node,GridNodePtr node_1)
{
    //id=0,计算score,放入开集
    if (node_1->id==0)
    {
        //从父节点更新
        node_1->id=1;
        node_1->cameFrom=node;
        node_1->gScore=node->gScore+1;
        node_1->fScore=fabs(node_1->coord(0)-target->coord(0))+
                       fabs(node_1->coord(1)-target->coord(1))+
                       fabs(node_1->coord(2)-target->coord(2));
        openSet.insert({node_1->gScore+node_1->fScore,node_1});
        return;
    }
    //id=1,计算score与现有score比较
    if (node_1->id==1)
    {
        //更新开集
        if (node_1->gScore>node->gScore+1)
        {
            node_1->cameFrom=node;
            node_1->gScore=node->gScore+1;
            openSet.insert({node_1->gScore+node_1->fScore,node_1});
        }
        return;
    }
    //id=-1，pass
    if (node_1->id==-1)
        return;
}

void gridPathFinder::showGridPathFun()
{
    showGridPath.push_back(target->coord);
    GridNodePtr node=target;

    while(node->cameFrom!=NULL)
    {
        showGridPath.push_back(node->cameFrom->coord);
        node=node->cameFrom;
    }

}