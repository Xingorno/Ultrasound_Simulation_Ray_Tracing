#ifndef OBJLOADER_H
#define OBJLOADER_H

#include "Bullet3Common/b3MinMax.h"

// load_mesh_from_obj
#include "wavefront/tiny_obj_loader.h"
#include "Bullet3Common/b3AlignedObjectArray.h"
#include <string>
#include <vector>

// btgCreateGraphicsShapeFromWavefrontObj
#include "btBulletDynamicsCommon.h"
#include "OpenGL/GLInstanceGraphicsShape.h"

GLInstanceGraphicsShape* btgCreateGraphicsShapeFromWavefrontObj(std::vector<tinyobj::shape_t>& shapes, bool flatShading=false)
{

    b3AlignedObjectArray<GLInstanceVertex>* vertices = new b3AlignedObjectArray<GLInstanceVertex>;
    {
        b3AlignedObjectArray<int>* indicesPtr = new b3AlignedObjectArray<int>;

        for (int s = 0; s<static_cast<int>(shapes.size()); s++)
        {
            tinyobj::shape_t& shape = shapes[s];
            int faceCount = shape.mesh.indices.size();

            for (int f=0;f<faceCount;f+=3)
            {
                btVector3 normal(0,1,0);
                int vtxBaseIndex = vertices->size();

                GLInstanceVertex vtx0;
                vtx0.xyzw[0] = shape.mesh.positions[shape.mesh.indices[f]*3+0];
                vtx0.xyzw[1] = shape.mesh.positions[shape.mesh.indices[f]*3+1];
                vtx0.xyzw[2] = shape.mesh.positions[shape.mesh.indices[f]*3+2];
                vtx0.xyzw[3] = 0.f;

                if (shape.mesh.texcoords.size())
                {
                    vtx0.uv[0] = shape.mesh.texcoords[shape.mesh.indices[f]*2+0];
                    vtx0.uv[1] = shape.mesh.texcoords[shape.mesh.indices[f]*2+1];
                }
                else
                {
                    vtx0.uv[0] = 0.5;
                    vtx0.uv[1] = 0.5;
                }

                GLInstanceVertex vtx1;
                vtx1.xyzw[0] = shape.mesh.positions[shape.mesh.indices[f+1]*3+0];
                vtx1.xyzw[1] = shape.mesh.positions[shape.mesh.indices[f+1]*3+1];
                vtx1.xyzw[2] = shape.mesh.positions[shape.mesh.indices[f+1]*3+2];
                vtx1.xyzw[3]= 0.f;

                if (shape.mesh.texcoords.size())
                {
                    vtx1.uv[0] = shape.mesh.texcoords[shape.mesh.indices[f+1]*2+0];
                    vtx1.uv[1] = shape.mesh.texcoords[shape.mesh.indices[f+1]*2+1];
                }
                else
                {
                    vtx1.uv[0] = 0.5f;
                    vtx1.uv[1] = 0.5f;
                }

                GLInstanceVertex vtx2;
                vtx2.xyzw[0] = shape.mesh.positions[shape.mesh.indices[f+2]*3+0];
                vtx2.xyzw[1] = shape.mesh.positions[shape.mesh.indices[f+2]*3+1];
                vtx2.xyzw[2] = shape.mesh.positions[shape.mesh.indices[f+2]*3+2];
                vtx2.xyzw[3] = 0.f;
                if (shape.mesh.texcoords.size())
                {
                    vtx2.uv[0] = shape.mesh.texcoords[shape.mesh.indices[f+2]*2+0];
                    vtx2.uv[1] = shape.mesh.texcoords[shape.mesh.indices[f+2]*2+1];
                } else
                {
                    vtx2.uv[0] = 0.5;
                    vtx2.uv[1] = 0.5;
                }

                btVector3 v0(vtx0.xyzw[0],vtx0.xyzw[1],vtx0.xyzw[2]);
                btVector3 v1(vtx1.xyzw[0],vtx1.xyzw[1],vtx1.xyzw[2]);
                btVector3 v2(vtx2.xyzw[0],vtx2.xyzw[1],vtx2.xyzw[2]);

                unsigned int maxIndex = 0;
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f]*3+0);
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f]*3+1);
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f]*3+2);
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f+1]*3+0);
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f+1]*3+1);
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f+1]*3+2);
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f+2]*3+0);
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f+2]*3+1);
                maxIndex = b3Max(maxIndex,shape.mesh.indices[f+2]*3+2);
                bool hasNormals = (shape.mesh.normals.size() && maxIndex<shape.mesh.normals.size() );

                if (flatShading || !hasNormals)
                {
                    normal = (v1-v0).cross(v2-v0);
                    btScalar len2 = normal.length2();
                    //skip degenerate triangles
                    if (len2 > SIMD_EPSILON)
                    {
                        normal.normalize();
                    }
                    else
                    {
                        normal.setValue(0,0,0);
                    }
                    vtx0.normal[0] = normal[0];
                    vtx0.normal[1] = normal[1];
                    vtx0.normal[2] = normal[2];
                    vtx1.normal[0] = normal[0];
                    vtx1.normal[1] = normal[1];
                    vtx1.normal[2] = normal[2];
                    vtx2.normal[0] = normal[0];
                    vtx2.normal[1] = normal[1];
                    vtx2.normal[2] = normal[2];
                }
                else
                {
                    vtx0.normal[0] = shape.mesh.normals[shape.mesh.indices[f]*3+0];
                    vtx0.normal[1] = shape.mesh.normals[shape.mesh.indices[f]*3+1];
                    vtx0.normal[2] = shape.mesh.normals[shape.mesh.indices[f]*3+2]; //shape.mesh.indices[f+1]*3+0
                    vtx1.normal[0] = shape.mesh.normals[shape.mesh.indices[f+1]*3+0];
                    vtx1.normal[1] = shape.mesh.normals[shape.mesh.indices[f+1]*3+1];
                    vtx1.normal[2] = shape.mesh.normals[shape.mesh.indices[f+1]*3+2];
                    vtx2.normal[0] = shape.mesh.normals[shape.mesh.indices[f+2]*3+0];
                    vtx2.normal[1] = shape.mesh.normals[shape.mesh.indices[f+2]*3+1];
                    vtx2.normal[2] = shape.mesh.normals[shape.mesh.indices[f+2]*3+2];
                }
                vertices->push_back(vtx0);
                vertices->push_back(vtx1);
                vertices->push_back(vtx2);
                indicesPtr->push_back(vtxBaseIndex);
                indicesPtr->push_back(vtxBaseIndex+1);
                indicesPtr->push_back(vtxBaseIndex+2);
            }
        }


        GLInstanceGraphicsShape* gfxShape = new GLInstanceGraphicsShape;
        gfxShape->m_vertices = vertices;
        gfxShape->m_numvertices = vertices->size();
        gfxShape->m_indices = indicesPtr;
        gfxShape->m_numIndices = indicesPtr->size();
        for (int i=0;i<4;i++)
            gfxShape->m_scaling[i] = 1;//bake the scaling into the vertices
        return gfxShape;
    }
}

GLInstanceGraphicsShape* load_mesh_from_obj(const std::string & relativeFileName, const std::string & materialPrefixPath)
{
    std::vector<tinyobj::shape_t> shapes;
    std::string err = tinyobj::LoadObj(shapes, relativeFileName.c_str(), materialPrefixPath.c_str());

    GLInstanceGraphicsShape* gfxShape = btgCreateGraphicsShapeFromWavefrontObj(shapes);
    return gfxShape;
}

#endif // OBJLOADER_H
