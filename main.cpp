/******************************************************************************
 * Ray tracing skeleton program.
 *
 * Happy hacking! - eric
 *****************************************************************************/

// required on macOS to compile the stb_image_write library without warnings
#ifdef __clang__
# pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

// include files
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <glm/gtc/type_ptr.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct Ray
{
    // ray parameters
    glm::vec3 origin;
    glm::vec3 direction;
};

struct Material
{
    // material components (vec3 variables are in RGB format)
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;
    float shininess;
};

struct SceneObject
{
    // the material for this scene object (with ambient, diffuse, specular, and shininess components)
    Material material;

    virtual ~SceneObject()
    {
        // the destructor needs to be declared as virtual to properly deallocate objects on exit
    }

    /**
     * @brief       Template method for calculating the intersection of this object with the provided ray
     * @param[in]   incomingRay            Ray that will be checked for intersection with this object
     * @param[out]  outIntersectionPoint   Point of intersection (in case there is an intersection)
     * @param[out]  outIntersectionNormal  Normal vector at the point of intersection (in case there is an intersection)
     * @return      If there is an intersection, returns the distance from the ray origin to the intersection point. Otherwise, returns a negative number.
     */
    virtual float Intersect(const Ray& incomingRay,
                            glm::vec3& outIntersectionPoint,
                            glm::vec3& outIntersectionNormal)
    {
        // a negative number is returned to indicate no intersection
        // (don't modify this here -- implement a subclass of SceneObject to customize this)
        return -INFINITY;
    }
};

// Subclass of SceneObject representing a Sphere scene object
struct Sphere : public SceneObject
{
    // sphere parameters
    glm::vec3 center;
    float radius;

    /**
     * @brief       Ray-sphere intersection
     * @param[in]   incomingRay            Ray that will be checked for intersection with this object
     * @param[out]  outIntersectionPoint   Point of intersection (in case there is an intersection)
     * @param[out]  outIntersectionNormal  Normal vector at the point of intersection (in case there is an intersection)
     * @return      If there is an intersection, returns the distance from the ray origin to the intersection point. Otherwise, returns a negative number.
     */
    virtual float Intersect(const Ray& incomingRay,
                            glm::vec3& outIntersectionPoint,
                            glm::vec3& outIntersectionNormal)
    {
         float t = -INFINITY;
        
        glm::vec3 l = center - incomingRay.origin;
        float tca = glm::dot(l, incomingRay.direction);
        if (tca < 0) return t;

        float d2 = glm::dot(l, l) - tca * tca;
        if (d2 > radius * radius) return t;

        float thc = sqrt(radius * radius - d2);
        float t0 = tca - thc;
        float t1 = tca + thc;

        if (t0 > t1) std::swap(t0, t1);
        if (t0 < 0) {
            t0 = t1;
            if (t0 < 0) return t;
        }

        outIntersectionPoint = incomingRay.origin + t0 * incomingRay.direction;
        outIntersectionNormal = glm::normalize(outIntersectionPoint - center);
        t = t0;

        return t;
    }
};

// Subclass of SceneObject representing a Triangle scene object
struct Triangle : public SceneObject
{
    // triangle vertices
    glm::vec3 A;
    glm::vec3 B;
    glm::vec3 C;

    /**
     * @brief       Ray-triangle intersection
     * @param[in]   incomingRay            Ray that will be checked for intersection with this object
     * @param[out]  outIntersectionPoint   Point of intersection (in case there is an intersection)
     * @param[out]  outIntersectionNormal  Normal vector at the point of intersection (in case there is an intersection)
     * @return      If there is an intersection, returns the distance from the ray origin to the intersection point. Otherwise, returns a negative number.
     */
    virtual float Intersect(const Ray& incomingRay,
                            glm::vec3& outIntersectionPoint,
                            glm::vec3& outIntersectionNormal)
    {
        // This is an implementation of the Möller-Trumbore intersection algorithm,
        // as shown in the slides, using code adapted from the original paper at
        // https://www.graphics.cornell.edu/pubs/1997/MT97.pdf.

            // calculate two triangle edges and the determinant
        float t = -INFINITY;
        glm::vec3 E1 = B - A;
        glm::vec3 E2 = C - A;
        glm::vec3 H = glm::cross(incomingRay.direction, E2);
        float det = glm::dot(H, E1);

        // if determinant is less than 0, the triangle is facing backwards;
        // if determinant is close to 0, the ray is parallel to the triangle
        // (either way, there is no intersection)
        if (det < glm::epsilon<float>())
            return t;  // return -INFINITY, i.e., no intersection

        // compute u; if out of range, the ray is not intersecting the triangle
        glm::vec3 S = incomingRay.origin - A;
        float u = glm::dot(H, S) / det;
        if (u < 0.0f || u > 1.0f)
            return t;  // return -INFINITY, i.e., no intersection

        // compute v; if out of range, the ray is not intersecting the triangle
        glm::vec3 Q = glm::cross(S, E1);
        float v = glm::dot(Q, incomingRay.direction) / det;
        if (v < 0.0f || u + v > 1.0f)
            return t;  // return -INFINITY, i.e., no intersection

        // compute the distance t from the origin of the ray to the intersection point
        t = glm::dot(Q, E2) / det;

        // calculate the intersection point and normal vector
        outIntersectionPoint = incomingRay.origin + t * incomingRay.direction;
        outIntersectionNormal = glm::normalize(glm::cross(E1, E2));

        return t;

    }
};

struct Camera
{
    // physical camera parameters
    glm::vec3 position;
    glm::vec3 lookTarget;
    glm::vec3 globalUp;
    float fovY;
    float focalLength;

    // final image parameters
    int imageWidth;
    int imageHeight;
};

struct Light
{
    // light position (w = 1 if point light, w = 0 if directional light)
    glm::vec4 position;

    // light components (in RGB)
    glm::vec3 ambient;
    glm::vec3 diffuse;
    glm::vec3 specular;

    // attenuation factors
    float constant;
    float linear;
    float quadratic;
};

struct IntersectionInfo
{
    Ray incomingRay;               // ray used to calculate the intersection
    float t;                       // distance from the ray's origin to the point of intersection
    SceneObject* object;           // object the ray intersected with (if == nullptr, there was no intersection)
    glm::vec3 intersectionPoint;   // point where the intersection occured
    glm::vec3 intersectionNormal;  // normal vector at the point of intersection
};

struct Scene
{
    std::vector<SceneObject*> objects;  // list of all objects in the scene
    std::vector<Light> lights;          // list of all lights in the scene
};

struct Image
{
    std::vector<unsigned char> data;  // image data
    int width;                        // image width
    int height;                       // image height

    /**
     * @brief      Constructor
     * @param[in]  w  Width
     * @param[in]  h  Height
     */
    Image(const int& w, const int& h)
        : width(w)
        , height(h)
    {
        data.resize(w * h * 3, 0);
    }

    /**
     * @brief      Converts the provided color value from [0, 1] to [0, 255]
     * @param[in]  c  Color value in [0, 1] range
     * @return     Color value in [0, 255] range
     */
    unsigned char ToChar(float c)
    {
        c = glm::clamp(c, 0.0f, 1.0f);
        return static_cast<unsigned char>(c * 255);
    }

    /**
     * @brief      Sets the color at the specified pixel location
     * @param[in]  x      X-coordinate of the pixel
     * @param[in]  y      Y-coordinate of the pixel
     * @param[in]  color  Pixel color
     */
    void SetColor(const int& x, const int& y, const glm::vec3& color)
    {
        int index = (y * width + x) * 3;
        data[index] = ToChar(color.r);
        data[index + 1] = ToChar(color.g);
        data[index + 2] = ToChar(color.b);
    }
};

/**
 * @brief      Gets the ray that goes from the camera's position to the specified pixel at (x, y)
 * @param[in]  camera  Camera data
 * @param[in]  x       X-coordinate of the pixel
 * @param[in]  y       Y-coordinate of the pixel
 * @return     Ray that passes through the pixel at (x, y)
 */
Ray GetRayThruPixel(const Camera &camera, const int& pixelX, const int& pixelY)
{
    // This function is implemented exactly as shown in the slides,
    // to give you a rough idea how you would implement the other
    // functions using the glm library.

    // get the width and height of the viewport
    float aspect = (float) camera.imageWidth / camera.imageHeight;
    float hViewport = 2 * camera.focalLength * glm::tan(glm::radians(camera.fovY) / 2);
    float wViewport = aspect * hViewport;

    // calculate the camera's basis vectors l, u, and v
    glm::vec3 l = glm::normalize(camera.lookTarget - camera.position);
    glm::vec3 u = glm::normalize(glm::cross(l, camera.globalUp));
    glm::vec3 v = glm::normalize(glm::cross(u, l));

    // calculate the lower-left corner point L of the viewport
    glm::vec3 L = camera.position + l * camera.focalLength - u * (wViewport / 2) - v * (hViewport / 2);

    // calculate the ray that goes through the center of the pixel (pixelX, pixelY)
    float s = (pixelX + 0.5f) / camera.imageWidth * wViewport;
    float t = (pixelY + 0.5f) / camera.imageHeight * hViewport;
    glm::vec3 P = L + u * s + v * t;

    // finally, construct the ray
    Ray ray;
    ray.origin = camera.position;
    ray.direction = glm::normalize(P - ray.origin);
    return ray;
}

/**
 * @brief      Cast a ray into the scene
 * @param[in]  ray    Ray to cast to the scene
 * @param[in]  scene  Scene object
 * @return     An IntersectionInfo object that will contain the results of the raycast
 */
IntersectionInfo RayCast(const Ray& ray, const Scene& scene)
{
    IntersectionInfo intersect;
    intersect.incomingRay = ray;
    intersect.t = INFINITY;
    intersect.object = nullptr;

    // this naive algorithm simply checks, for all scene objects,
    // whether a given scene object intersects with the ray
    for (size_t i = 0; i < scene.objects.size(); i++)
    {
        glm::vec3 point, normal;
        float t = scene.objects[i]->Intersect(ray, point, normal);

        // if there is an intersection, and it is closer than a previous scene object's intersection...
        if (t >= 0.0f && t < intersect.t)
        {
            // ... store that intersection into the IntersectionInfo object that we will return
            intersect.t = t;
            intersect.object = scene.objects[i];
            intersect.intersectionPoint = point;
            intersect.intersectionNormal = normal;
        }
    }

    return intersect;
}

/**
 * @brief      Perform a ray-trace to the scene
 * @param[in]  ray        Ray to trace
 * @param[in]  scene      Scene data
 * @param[in]  cameraPos  Original camera position (needed for computing specular highlights)
 * @param[in]  depth      Current reflection depth of this trace (use this to track your recursion)
 * @return     Resulting color after the ray bounced around the scene
 */
glm::vec3 RayTrace(const Ray& ray, const Scene& scene, glm::vec3& cameraPos, int depth)
{
    glm::vec3 color(0.0f);

    // cast a ray into the scene
    IntersectionInfo intersect = RayCast(ray, scene);

    // if there was an intersection...
    if (intersect.object != nullptr)
    {
        Material material = intersect.object->material;
        glm::vec3 ambient = material.ambient;
        glm::vec3 diffuse = material.diffuse;
        glm::vec3 specular = material.specular;
        float shininess = material.shininess;
        glm::vec3 normal = intersect.intersectionNormal;
        glm::vec3 intersectionPoint = intersect.intersectionPoint;
        glm::vec3 viewDir = glm::normalize(cameraPos - intersectionPoint);

        // calculate the color of each light
        for (Light light : scene.lights)
        {
            glm::vec3 lightDir;
            float attenuation = 1.0f;

            // calculate the light direction and attenuation based on whether the light is a point or directional light
            if (light.position.w == 1.0f)
            {
                lightDir = glm::normalize(glm::vec3(light.position) - intersectionPoint);
                float distance = glm::distance(intersectionPoint, glm::vec3(light.position));
                attenuation = 1.0f / (light.constant + light.linear * distance + light.quadratic * (distance * distance));
            }
            else
            {
                lightDir = glm::normalize(-glm::vec3(light.position));
            }

            // calculate the ambient, diffuse, and specular components for this light
            glm::vec3 ambientComp = ambient * light.ambient;
            glm::vec3 diffuseComp = diffuse * light.diffuse * glm::max(glm::dot(normal, lightDir), 0.0f);
            glm::vec3 reflectDir = glm::normalize(glm::reflect(-lightDir, normal));
            glm::vec3 specularComp = specular * light.specular * glm::pow(glm::max(glm::dot(viewDir, reflectDir), 0.0f), shininess);

            // add this light's contribution to the color
            color += attenuation * (ambientComp + diffuseComp + specularComp);
        }
    }

    return color;
}



/**
 * @brief       Read a .test file into our Scene structures
 * @param[in]   filename  Filename of the scene to render
 * @param[out]  camera    Camera structure
 * @param[out]  maxDepth  Maximum recursion depth
 * @param[out]  scene     Scene structure to populate with objects
 * @return      true if the scene was successfully load
 */
bool LoadScene(char* filename, Camera& camera, int& maxDepth, Scene& scene)
{
    std::ifstream file(filename);
    if (! file)
        return false;

    // populate Camera parameters
    file >> camera.imageWidth;
    file >> camera.imageHeight;
    file >> camera.position.x;
    file >> camera.position.y;
    file >> camera.position.z;
    file >> camera.lookTarget.x;
    file >> camera.lookTarget.y;
    file >> camera.lookTarget.z;
    file >> camera.globalUp.x;
    file >> camera.globalUp.y;
    file >> camera.globalUp.z;
    file >> camera.fovY;
    file >> camera.focalLength;

    // read the max depth
    file >> maxDepth;

    // read the number of scene objects
    size_t n;
    file >> n;

    // for each object...
    for (size_t i = 0; i < n; i++)
    {
        // read the type of object
        std::string objectType;
        file >> objectType;

        SceneObject* object;
        if (objectType == "sphere")
        {
            // create a new sphere
            Sphere* sphere = new Sphere();
            file >> sphere->center.x;
            file >> sphere->center.y;
            file >> sphere->center.z;
            file >> sphere->radius;
            object = sphere;
        }
        else if (objectType == "tri")
        {
            // create a new triangle
            Triangle* triangle = new Triangle();
            file >> triangle->A.x;
            file >> triangle->A.y;
            file >> triangle->A.z;
            file >> triangle->B.x;
            file >> triangle->B.y;
            file >> triangle->B.z;
            file >> triangle->C.x;
            file >> triangle->C.y;
            file >> triangle->C.z;
            object = triangle;
        }
        else
        {
            std::cout << "Invalid object type " << objectType << "\n";
            return false;
        }

        // populate the Material properties
        file >> object->material.ambient.r;
        file >> object->material.ambient.g;
        file >> object->material.ambient.b;
        file >> object->material.diffuse.r;
        file >> object->material.diffuse.g;
        file >> object->material.diffuse.b;
        file >> object->material.specular.r;
        file >> object->material.specular.g;
        file >> object->material.specular.b;
        file >> object->material.shininess;

        scene.objects.push_back(object);
    }

    // read the number of lights
    file >> n;

    // for each light...
    for (size_t i = 0; i < n; i++)
    {
        Light light;

        // populate the Light properties
        file >> light.position.x;
        file >> light.position.y;
        file >> light.position.z;
        file >> light.position.w;
        file >> light.ambient.r;
        file >> light.ambient.g;
        file >> light.ambient.b;
        file >> light.diffuse.r;
        file >> light.diffuse.g;
        file >> light.diffuse.b;
        file >> light.specular.r;
        file >> light.specular.g;
        file >> light.specular.b;
        file >> light.constant;
        file >> light.linear;
        file >> light.quadratic;

        scene.lights.push_back(light);
    }

    return true;
}

/**
 * Main function
 */
int main(int argc, char* argv[])
{
    // load the scene
    Scene scene;
    Camera camera;
    int maxDepth;
    if (argc < 2)
    {
        std::cout << "Specify a .test file to render.\n";
        return 1;
    }
    bool success = LoadScene(argv[1], camera, maxDepth, scene);
    if (! success)
    {
        std::cout << "Error reading " << argv[1] << "\n";
        return 1;
    }

    // declare an image object
    Image image(camera.imageWidth, camera.imageHeight);

    // for each pixel in the image, trace a ray into the scene to determine the pixel's color
    for (int y = 0; y < image.height; y++)
    {
        for (int x = 0; x < image.width; x++)
        {
            Ray ray = GetRayThruPixel(camera, x, image.height - y - 1);

            glm::vec3 color = RayTrace(ray, scene, camera.position, maxDepth);
            image.SetColor(x, y, color);
        }

        // update user with the current row being drawn
        std::cout << "Row: "
                  << std::setfill(' ') << std::setw(4) << (y + 1)
                  << " / "
                  << std::setfill(' ') << std::setw(4) << image.height
                  << "\r" << std::flush;
    }
    std::cout << std::endl;

    // write the image as a PNG file
    stbi_write_png("test.png", image.width, image.height, 3, image.data.data(), 0);

    // cleanup
    for (size_t i = 0; i < scene.objects.size(); i++)
        delete scene.objects[i];

    return 0;
}
