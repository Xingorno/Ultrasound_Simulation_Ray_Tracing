#include "scene.h"

#include "objloader.h"
#include "ray.h"

#include <cmath>
#include <iostream>
#include <cassert>
#include <exception>

//template std::array<std::vector<ray_physics::segment>, 256> scene::cast_rays<256>();
// template std::array<std::vector<ray_physics::segment>, 1024> scene::cast_rays<1024>();
// template std::array<std::vector<ray_physics::segment>, 2048> scene::cast_rays<2048>();
template std::array<std::array<std::vector<ray_physics::segment>, 20>, 1024> scene::cast_rays<20, 1024>(transducer_ & transducer);
// template std::array<std::vector<ray_physics::segment>, 512> scene::cast_rays<512>();
scene::scene(const nlohmann::json & config, transducer_ & transducer) :
    transducer(transducer)
{
    try
    {
        parse_config(config);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error{ "Error while loading scene: " + std::string{ex.what()} };
    }

    create_empty_world();

    init();
}

scene::~scene()
{
    destroy_world();
}

void scene::init()
{
    for (auto & mesh : meshes)
    {
        const auto full_path = working_dir + mesh.filename;

        auto object = add_rigidbody_from_obj(full_path, mesh.deltas, scaling); // delta relates to the position of objects; scaling realtes to the scale of object in rendering scene

        object->setUserPointer(&mesh);
    }
}

// template<unsigned int ray_count>
// std::array<std::vector<ray_physics::segment>, ray_count> scene::cast_rays()
template<unsigned int sample_count, unsigned int ray_count>
std::array<std::array<std::vector<ray_physics::segment>, sample_count>, ray_count> scene::cast_rays(transducer_ & transducer)
{
    using namespace ray_physics;

    // std::array<std::vector<segment>, ray_count> segments;
    std::array<std::array<std::vector<segment>, sample_count>, ray_count> segments;

    unsigned int tests = 0;
    unsigned int total_collisions = 0;

    ///step the simulation
    if (m_dynamicsWorld)
    {
        const float ray_start_step { 2.0f }; // move transducer along fixed gap [mm]

        // for (auto & segments_vector : segments)
        // {
        //     segments_vector.reserve(std::pow(2, ray::max_depth));
        // }
        // allocate vecotr capacity
        for (auto & sample_vector : segments) //understanding: https://stackoverflow.com/questions/35490236/colon-and-auto-in-for-loop-c-need-some-help-understanding-the-syntax
        {
            for (auto & segments_vector : sample_vector )
            {
                segments_vector.reserve(ray::max_depth);
            }
        }
        size_t ray_i = 473;
        // size_t ray_i = 343;
        for (size_t ray_i = 0; ray_i < ray_count; ray_i++)
        {
            // auto & segments_vector = segments[ray_i];
            auto & samples_vector = segments[ray_i];
            for(size_t sample_i=0; sample_i < sample_count; sample_i++)
            {
                std::vector<ray> ray_stack;
                ray_stack.reserve(ray::max_depth-1);
                // std::array<ray, sample_count> samples;
                // Add first ray
                // {
                ray first_ray
                {
                    transducer.element(ray_i, sample_i).position,                          // from [mm]
                    transducer.element(ray_i, sample_i).direction,                         // initial direction
                    0,                                                           // depth
                    materials.at(starting_material),
                    nullptr,
                    initial_intensity,
                    transducer.frequency,
                    units::length::millimeter_t(0),                              // distance traveled [mm]
                    0                                                            // previous ray
                };

                ray_stack.push_back(first_ray);
                // for(size_t sample_i = 0; sample_i < sample_count; sample_i++)
                // {
                //     samples.at(sample_i) = first_ray;
                //     
                // }
                // }
                // indexMove++;
                while (ray_stack.size() > 0)
                // for (size_t sample_i = 0; sample_i < sample_count; sample_i++)
                {
               
                    // Pop a ray from the stack and check if it collides
                     // auto & ray_ = samples.at(sample_i);
                    auto & ray_ = ray_stack.at(ray_stack.size()-1);
                    // if (!ray_.null)
                    // {
                    float r_length = ray_physics::max_ray_length(ray_); //[mm]
                    auto to = ray_.from + enlarge(ray_.direction, r_length); //[mm]
            
                    btCollisionWorld::ClosestRayResultCallback closestResults(ray_.from + ray_.direction ,to); // adding ray_direction * 0.1 is for uncorrect ray tracing operation nearby the "from" point 
                    
                    m_dynamicsWorld->rayTest(ray_.from + ray_.direction,to,closestResults);

                    btCollisionWorld::ClosestRayResultCallback closestResults1(closestResults.m_hitPointWorld,to);
                    m_dynamicsWorld->rayTest(closestResults.m_hitPointWorld,to,closestResults1);
  
                    
                    tests++;

                    ray_stack.pop_back();
                    if (closestResults.hasHit())
                    {
                        if (ray_.depth < ray::max_depth)
                        {
                        // Substract ray intensity according to distance traveled
                        auto distance_before_hit = ray_.distance_traveled;
                        auto intensity_before_hit = ray_.intensity;
                        
                        ray_physics::travel(ray_, distance_in_mm(ray_.from, closestResults.m_hitPointWorld));

                        // Calculate refraction and reflection directions and intensities
                        //compute the distance between the closet point and all hit point
                        mesh* organ;
                        mesh* organ1;
                        if (!closestResults1.hasHit())
                        {
                            organ = static_cast<mesh*>(closestResults.m_collisionObject->getUserPointer());
                        }
                        else
                        {
                            float distance = closestResults.m_hitPointWorld.distance(closestResults1.m_hitPointWorld);
                            if ( distance < 1.2)
                            {   
                                btVector3 delta = closestResults.m_hitPointWorld - closestResults1.m_hitPointWorld;
                                float vector_product = ray_.direction.dot(delta);
                                
                                if (vector_product < 0)
                                {
                                    organ = static_cast<mesh*>(closestResults1.m_collisionObject->getUserPointer());
                                }
                                else
                                {
                                    organ = static_cast<mesh*>(closestResults.m_collisionObject->getUserPointer());
                                }
                                
                                // organ = static_cast<mesh*>(closestResults1.m_collisionObject->getUserPointer());
                                organ1= static_cast<mesh*>(closestResults.m_collisionObject->getUserPointer());
                            }
                            else
                            {
                                organ = static_cast<mesh*>(closestResults.m_collisionObject->getUserPointer());
                            }
                        }
                        
                        
                        
                        // const auto organ= static_cast<mesh*>(closestResults.m_collisionObject->getUserPointer());
                        assert(organ);

                        auto result = ray_physics::hit_boundary(ray_, closestResults.m_hitPointWorld, closestResults.m_hitNormalWorld, *organ);

                        // Register collision creating a segment from the beginning of the ray to the collision point
                        //segments_vector.emplace_back(segment{ray_.from, closestResults.m_hitPointWorld, ray_.direction, result.reflected_intensity, intensity_before_hit, ray_.media.attenuation, result.reflection.distance_traveled, ray_.media});
                        

                        // segments_vector.emplace_back(segment{ray_.from, closestResults.m_hitPointWorld, ray_.direction, result.reflected_intensity, intensity_before_hit, ray_.media.attenuation, distance_before_hit, ray_.media});
                        samples_vector[sample_i].emplace_back(segment{ray_.from, closestResults.m_hitPointWorld, ray_.direction, result.reflected_intensity, intensity_before_hit, ray_.media.attenuation, distance_before_hit, ray_.media});
                        // Spawn reflection and refraction rays
                       
                        if (result.refraction.intensity > ray::intensity_epsilon)
                        // if (result.refraction.intensity > 0.0001)
                        {   
                            // result.refraction.parent_collision = segments_vector.size()-1;
                            result.refraction.parent_collision = samples_vector[sample_i].size()-1;
                            ray_stack.push_back(result.refraction);
                            // samples.at(sample_i) = result.refraction;
                        }

                        // if (result.reflection.intensity > ray::intensity_epsilon)
                        // {
                        //     result.reflection.parent_collision = segments_vector.size()-1;
                        //     ray_stack.push_back(result.reflection);
                        // }
                        }
                    
                    }
                    else
                    {
                    // Ray did not reach another media, add a data point at its end.
                    // segments_vector.emplace_back(segment{ray_.from, to, ray_.direction, 0.0f, ray_.intensity, ray_.media.attenuation, ray_.distance_traveled, ray_.media});
                        auto to_1 = ray_.from + enlarge(ray_.direction, 150);
                        samples_vector[sample_i].emplace_back(segment{ray_.from, to_1, ray_.direction, 0.0f, ray_.intensity, ray_.media.attenuation, ray_.distance_traveled, ray_.media});
                    // samples.at(sample_i).null = true;
                    }
                }  
            }
        }
              
        for (auto & collision_vector : segments)
        {
            total_collisions += collision_vector.size();
        }
    }

    const float fps = 1.0f/(float( clock() - frame_start ) /  CLOCKS_PER_SEC);
    std::cout <<"fps: "<< fps << "; number of tests: " << tests << "; number of collisions: " << total_collisions << std::endl;
    frame_start = clock();

    return segments;
}

void scene::parse_config(const nlohmann::json & config)
{
    using namespace units::angle;

    working_dir = config.find("workingDirectory") != config.end() ? config.at("workingDirectory") : "";

    const auto & t_pos = config.at("transducerPosition");
    transducer_pos = {t_pos[0], t_pos[1], t_pos[2]};

    const auto & orig = config.at("origin");
    origin = {orig[0], orig[1], orig[2]};

    const auto & spac = config.at("spacing");
    spacing = {spac[0], spac[1], spac[2]};

    starting_material = config.at("startingMaterial");

    scaling = config.at("scaling");

    const auto & mats = config.at("materials");
    if (mats.is_array())
    {
        for (const auto & mat : mats)
        {
            materials[mat.at("name")] =
                {
                    mat.at("impedance"),
                    mat.at("attenuation"),
                    mat.at("mu0"),
                    mat.at("mu1"),
                    mat.at("sigma"),
                    mat.at("specularity"),
                    mat.at("roughness")
                };
        }
    }
    else
    {
        throw std::runtime_error("materials must be an array");
    }

    const auto & meshes_ = config.at("meshes");
    if (meshes_.is_array())
    {
        for (const auto & mesh_ : meshes_)
        {
            const auto & deltas { mesh_.at("deltas") };
            meshes.emplace_back(mesh{
                        mesh_.at("file"),
                        mesh_.at("rigid"),
                        mesh_.at("vascular"),
                        {deltas[0], deltas[1], deltas[2]},
                        mesh_.at("outsideNormals"),
                        materials.at(mesh_.at("material")),
                        materials.at(mesh_.at("outsideMaterial"))});
        }
    }
    else
    {
        throw std::runtime_error("meshes must be an array");
    }
}

void scene::create_empty_world()
{
    //collision configuration contains default setup for memory, collision setup
    m_collisionConfiguration = std::make_unique<btDefaultCollisionConfiguration>();

    //use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
    m_dispatcher = std::make_unique<btCollisionDispatcher>(m_collisionConfiguration.get());

    m_broadphase = std::make_unique<btDbvtBroadphase>();
    //the default constraint solver. For parallel processing you can use a different solver (see Extras/BulletMultiThreaded)
    m_solver = std::make_unique<btSequentialImpulseConstraintSolver>();

    m_dynamicsWorld = std::make_unique<btDiscreteDynamicsWorld>(m_dispatcher.get(),m_broadphase.get(),m_solver.get(),m_collisionConfiguration.get());

    m_dynamicsWorld->setGravity(btVector3(0,-10,0));
}

void scene::destroy_world()
{
    //delete collision shapes
    for (int j = 0; j < m_collisionShapes.size(); j++)
    {
        btCollisionShape* shape = m_collisionShapes[j];
        delete shape;
    }
    m_collisionShapes.clear();

    m_dynamicsWorld.reset();
    m_solver.reset();
    m_broadphase.reset();
    m_dispatcher.reset();
    m_collisionConfiguration.reset();
}

units::length::millimeter_t scene::distance_in_mm(const btVector3 & v1, const btVector3 & v2) const
{
    using namespace std;

    auto x_dist = abs(v1.getX() - v2.getX()) * spacing[0];
    auto y_dist = abs(v1.getY() - v2.getY()) * spacing[1];
    auto z_dist = abs(v1.getZ() - v2.getZ()) * spacing[2];

    //mm
    return units::length::millimeter_t(sqrt(pow(x_dist,2) + pow(y_dist,2) + pow(z_dist,2)));
    //return units::length::millimeter_t(sqrt(pow(x_dist,2) + pow(y_dist,2) + pow(z_dist,2)) * 10);
}

btVector3 scene::enlarge(const btVector3 & versor, float mm) const
{
    // assert(versor.length2() < 1.1f);
   
    return mm * btVector3 ( spacing[0] * versor.getX(),
                                   spacing[1] * versor.getY(),
                                   spacing[2] * versor.getZ() );

    // return mm/100.0f * btVector3 ( spacing[0] * versor.getX(),
    //                                spacing[1] * versor.getY(),
    //                                spacing[2] * versor.getZ() );
}

btRigidBody * scene::add_rigidbody_from_obj(const std::string & fileName, std::array<float, 3> deltas, float scaling)
{
    GLInstanceGraphicsShape* glmesh = load_mesh_from_obj(fileName, "");
    // printf("[INFO] Obj loaded: Extracted %d verticed from obj file [%s]\n", glmesh->m_numvertices, fileName);

    const GLInstanceVertex& v = glmesh->m_vertices->at(0);
    btTriangleIndexVertexArray* tiva = new btTriangleIndexVertexArray(glmesh->m_numIndices / 3, &glmesh->m_indices->at(0), 3* sizeof(int),
                                                                      glmesh->m_numvertices, (btScalar*)(&(v.xyzw[0])), sizeof(GLInstanceVertex));

    btBvhTriangleMeshShape* shape = new btBvhTriangleMeshShape(tiva, true);

    m_collisionShapes.push_back(shape);

    float _scaling[4] = {scaling,scaling,scaling,1};

    btVector3 localScaling(_scaling[0],_scaling[1],_scaling[2]);
    shape->setLocalScaling(localScaling);

    btTransform startTransform;
    startTransform.setIdentity();

    //std::array<float, 3> origin { -18, -22, -5 }; // origin for organs scene
    float pos[4] = {deltas[0]*_scaling[0]*_scaling[0],deltas[1]*_scaling[1]*_scaling[1],deltas[2]*_scaling[2]*_scaling[2],0};
    btVector3 position(pos[0] + origin[0], pos[1] + origin[1], pos[2] + origin[2]);
    startTransform.setOrigin(position);

    btScalar mass(0.f);
    btVector3 localInertia(0, 0, 0);
    auto myMotionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo cInfo(mass, myMotionState, shape, localInertia);
    auto body = new btRigidBody(cInfo);
    body->setUserIndex(-1);
    m_dynamicsWorld->addRigidBody(body);
    return body;
}

void scene::step(float delta_time)
{   // TODO: what does it mean?
    m_dynamicsWorld->stepSimulation(delta_time);
}


units::length::millimeter_t scene::distance(const btVector3 & from, const btVector3 & to) const
{
    return units::length::millimeter_t(from.distance(to));//[mm]
    //return units::length::millimeter_t(from.distance(to)*10.0f);
}
