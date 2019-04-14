
#include <stdio.h>
#include <sys/time.h>
#include <fstream>
#include <iostream>
#include <locale>
#include <map>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include <png.h>

#define TINYPLY_IMPLEMENTATION
#include "tinyply.h"
using namespace tinyply;

#include "color.h"
#include "la.h"
#include "light.h"
#include "material.h"
#include "scene.h"
#include "shape.h"
/*#include "ray.h"
#include "camera.h"
#include "sensor.h"

#include "tiny_obj_loader.cc"
*/

/*
    Transform state
*/

extern std::shared_ptr<Scene> scene;
extern Transform cameraTransform;

struct Params {
    std::unordered_map<std::string, std::vector<std::string>> strings;
    std::unordered_map<std::string, std::vector<float>> reals;
    std::unordered_map<std::string, std::vector<int64_t>> ints;
    std::unordered_map<std::string, std::vector<bool>> bools;

    bool hasString(std::string const& key) const {
        return strings.count(key) > 0;
    }

    std::string getString(std::string const& key) const {
        auto i = strings.find(key);
        if (i != strings.end())
            return i->second.at(0);
        else
            return "";
    }

    bool hasFloat(std::string const& key) const { return reals.count(key) > 0; }

    float getFloat(std::string const& key) const {
        auto i = reals.find(key);
        if (i != reals.end() && i->second.size() == 1)
            return i->second[0];
        else
            return 0.f;
    }

    bool hasRGB(std::string const& key) const { return reals.count(key) > 0; }

    rgba getRGB(std::string const& key) const {
        auto i = reals.find(key);
        if (i != reals.end() && i->second.size() == 3)
            return rgba(i->second[0], i->second[1], i->second[2], 1);
        else
            return rgba(0, 1, 0, 0);
    }
};

struct Object {
    std::vector<Instance> instances;
    std::vector<LightInstance> lights;
};

namespace PBRTState {

enum ActiveTransform {
    StartTime = 1,
    EndTime = 2,
    All = 3,
};

struct TransformState {
    float t0, t1;
    Transform start;
    Transform end;
    ActiveTransform active;
};

std::stack<TransformState> transforms;
std::unordered_map<std::string, TransformState> namedTransforms;

Transform cameraTransform;
std::vector<Shape const*> shapes;
// Shape const* worldShapes;

std::stack<std::pair<std::string, Object>> objectStack;
std::unordered_map<std::string, Object> namedObjects;

std::stack<std::shared_ptr<const Material>> materialStack;
std::unordered_map<std::string, std::shared_ptr<const Material>> namedMaterials;

std::unordered_map<std::string, std::shared_ptr<const Texture>> namedTextures;

std::unordered_map<std::string, std::shared_ptr<const Texture>> cachedTextures;
std::unordered_map<std::string, std::shared_ptr<const Shape>> cachedShapes;

void tSet(const Transform& t) {
    auto& top = transforms.top();
    if (top.active & ActiveTransform::StartTime) top.start = t;
    if (top.active & ActiveTransform::EndTime) top.end = t;
}

void tMul(const Transform& t) {
    auto& top = transforms.top();
    if (top.active & ActiveTransform::StartTime) top.start = t * top.start;
    if (top.active & ActiveTransform::EndTime) top.end = t * top.end;
}

struct SensorState {
    // Camera
    float shutteropen;
    float shutterclose;

    // Film
    uint32_t xresolution;
    uint32_t yresolution;
    float cropwindow[4];
    std::string filename;

    // Filter
    float xwidth;
    float ywidth;
};

SensorState sensor{

};

}  // namespace PBRTState

extern "C" {

std::string pbrtPath;

void pbrtInit(std::string path) {
    pbrtPath = path;
    PBRTState::transforms = {};
    PBRTState::namedTransforms.clear();
    PBRTState::transforms.push(PBRTState::TransformState{
        0, 1, Transform::Identity(), Transform::Identity(),
        PBRTState::ActiveTransform::All});
    PBRTState::shapes.clear();
    PBRTState::objectStack = {};
    PBRTState::namedObjects.clear();

    PBRTState::namedMaterials.clear();
    PBRTState::materialStack = {};
    PBRTState::materialStack.push(std::make_shared<DiffuseMaterial>(
        std::make_shared<ConstantTexture>(rgba(0, 0, 0, 0))));
    PBRTState::namedTextures.clear();

    PBRTState::cachedTextures.clear();
}

void pbrtCleanup() {
    PBRTState::transforms = {};
    PBRTState::namedTransforms.clear();
    PBRTState::shapes.clear();
    PBRTState::objectStack = {};
    PBRTState::namedObjects.clear();
    PBRTState::materialStack = {};
    PBRTState::namedTextures.clear();

    PBRTState::cachedTextures.clear();
}

// Transformations

void pbrtIdentity() { PBRTState::tSet(Transform::Identity()); }

void pbrtTranslate(float dx, float dy, float dz) {
    PBRTState::tMul(Transform::Translate(dx, dy, dz));
}

void pbrtTransform(float t[16]) {
    PBRTState::tSet(Transform(
        Point(t[0], t[4], t[8], t[12]), Point(t[1], t[5], t[9], t[13]),
        Point(t[2], t[6], t[10], t[14]), Point(t[3], t[7], t[11], t[15])));
}

void pbrtConcatTransform(float t[16]) {
    PBRTState::tMul(Transform(
        Point(t[0], t[4], t[8], t[12]), Point(t[1], t[5], t[9], t[13]),
        Point(t[2], t[6], t[10], t[14]), Point(t[3], t[7], t[11], t[15])));
}

void pbrtRotate(float angle, float x, float y, float z) {
    PBRTState::tMul(Transform::Rotate(V3(x, y, z), (angle / 180.0f) * M_PI));
}

void pbrtScale(float sx, float sy, float sz) {
    PBRTState::tMul(Transform::Scale(sx, sy, sz));
}

void pbrtLookAt(float ex, float ey, float ez, float lx, float ly, float lz,
                float ux, float uy, float uz) {
    PBRTState::tMul(
        Transform::LookAt(P3(ex, ey, ez), P3(lx, ly, lz), V3(ux, uy, uz)));
}

void pbrtCoordinateSystem(const char* name) {
    PBRTState::namedTransforms[std::string(name)] = PBRTState::transforms.top();
}

void pbrtCoordSysTransform(const char* name) {
    auto t = PBRTState::namedTransforms.find(std::string(name));
    if (t != PBRTState::namedTransforms.end()) {
        PBRTState::transforms.top() = t->second;
    }
}

void pbrtActiveTransformAll() {
    PBRTState::transforms.top().active = PBRTState::ActiveTransform::All;
}

void pbrtActiveTransformStartTime() {
    PBRTState::transforms.top().active = PBRTState::ActiveTransform::StartTime;
}

void pbrtActiveTransformEndTime() {
    PBRTState::transforms.top().active = PBRTState::ActiveTransform::EndTime;
}

void pbrtTransformTimes(float start, float end) {
    PBRTState::transforms.top().t0 = start;
    PBRTState::transforms.top().t1 = end;
}

void pbrtTransformBegin(void) {
    PBRTState::transforms.push(PBRTState::transforms.top());
}

void pbrtTransformEnd(void) { PBRTState::transforms.pop(); }
}

void pbrtCamera(std::string const& name, Params const& params) {
    PBRTState::cameraTransform = PBRTState::transforms.top().start;
    cameraTransform = PBRTState::transforms.top().start;
}

/*void pbrtSampler(std::string const& name, Params const& params) {
}

void pbrtFilm(std::string const& name, Params const& params) {
}

void pbrtPixelFilter(std::string const& name, Params const& params) {
}

void pbrtIntegrator(std::string const& name, Params const& params) {
}

void pbrtAccelerator(std::string const& name, Params const& params) {
}

void pbrtMakeNamedMedia(std::string const& name, Params const& params) {
}

void pbrtMediumInterface(std:string const& in, std::string const& out) {
}
*/

void pbrtAttributeBegin() {
    pbrtTransformBegin();
    PBRTState::materialStack.push(PBRTState::materialStack.top());
}

void pbrtAttributeEnd() {
    pbrtTransformEnd();
    PBRTState::materialStack.pop();
}

std::shared_ptr<const Shape> readply(std::string const& filepath) {
    // std::cout << pbrtPath << filepath << std::endl;

    PlyFile file;
    std::shared_ptr<PlyData> vertices, normals, uvs, faces;
    {
        std::ifstream ss(filepath, std::ios::binary);
        if (ss.fail()) throw std::runtime_error("failed to open " + filepath);

        file.parse_header(ss);

        try {
            vertices =
                file.request_properties_from_element("vertex", {"x", "y", "z"});
        } catch (const std::exception& e) {
            std::cout << "tinyply exception: " << e.what() << std::endl;
        }

        try {
            normals = file.request_properties_from_element("vertex",
                                                           {"nx", "ny", "nz"});
        } catch (const std::exception& e) { /*std::cout << "tinyply exception: "
                                               << e.what() << std::endl;*/
        }

        try {
            uvs = file.request_properties_from_element("vertex", {"u", "v"});
        } catch (const std::exception& e) { /*std::cout << "tinyply exception: "
                                               << e.what() << std::endl;*/
        }

        try {
            faces = file.request_properties_from_element("face",
                                                         {"vertex_indices"}, 3);
        } catch (const std::exception& e) {
            std::cout << "tinyply exception: " << e.what() << std::endl;
        }

        file.read(ss);
    }

    const size_t numVerticesBytes = vertices->buffer.size_bytes();
    std::vector<P3> vertexes(vertices->count);
    std::memcpy(vertexes.data(), vertices->buffer.get(), numVerticesBytes);

    std::vector<V3> n;
    if (normals) {
        const size_t numNormalsBytes = normals->buffer.size_bytes();
        n.resize(normals->count);
        std::memcpy(n.data(), normals->buffer.get(), numNormalsBytes);
    }

    std::vector<P2> u;
    if (uvs) {
        const size_t numUVBytes = uvs->buffer.size_bytes();
        u.resize(uvs->count);
        std::memcpy(u.data(), uvs->buffer.get(), numUVBytes);
    }

    const size_t numFacesBytes = faces->buffer.size_bytes();
    std::vector<Mesh::Triangle> triangles(faces->count);
    std::memcpy(triangles.data(), faces->buffer.get(), numFacesBytes);

    return std::make_shared<Mesh>(std::move(vertexes), std::move(n),
                                  std::move(u), std::move(triangles));
}

void pbrtShape(std::string const& name, Params const& params) {
    if (name == "plymesh") {
        auto filename = params.getString("string filename");
        auto cached = PBRTState::cachedShapes.find(pbrtPath + filename);
        std::shared_ptr<Shape const> shape;
        if (cached != PBRTState::cachedShapes.end()) {
            shape = cached->second;
        } else {
            shape = readply(pbrtPath + filename);
            PBRTState::cachedShapes[pbrtPath + filename] = shape;
        }

        PBRTState::objectStack.top().second.instances.push_back(
            Instance{shape, PBRTState::materialStack.top(),
                     PBRTState::transforms.top().start});
    } else {
        std::cout << "Unrecognized Shape: " << name << std::endl;
    }
}

void pbrtObjectBegin(std::string const& name) {
    PBRTState::objectStack.push({name, {}});
}

void pbrtObjectEnd() {
    if (PBRTState::objectStack.size() > 0) {
        auto object = PBRTState::objectStack.top();
        PBRTState::objectStack.pop();
        PBRTState::namedObjects[object.first] = object.second;
        std::cout << "Made object: " << object.first << " with "
                  << object.second.instances.size() << std::endl;
    }
}

void pbrtObjectInstance(std::string const& name) {
    auto const& object = PBRTState::namedObjects.at(name);

    // std::cout << "Instantiating " << name << " with " << shapes.size() <<
    // std::endl;
    // Instantiate by copying into current object stack with composed transform
    for (auto const& s : object.instances) {
        PBRTState::objectStack.top().second.instances.push_back(
            Instance{s.shape, s.material,
                     PBRTState::transforms.top().start * s.transform});
    }

    for (auto const& s : object.lights) {
        PBRTState::objectStack.top().second.lights.push_back(LightInstance{
            s.light, PBRTState::transforms.top().start * s.transform});
    }
}

std::shared_ptr<Texture const> makeTexture(std::string const& name,
                                           Params const& params,
                                           rgba const& defaultValue) {
    if (params.hasRGB(std::string("rgb ") + name)) {
        auto rgb = params.getRGB(std::string("rgb " + name));
        return std::make_shared<ConstantTexture const>(rgb);
    } else if (params.hasFloat(std::string("float ") + name)) {
        auto f = params.getFloat(std::string("float " + name));
        return std::make_shared<ConstantTexture const>(rgba(f, f, f, f));
    } else if (params.hasString(std::string("spectrum ") + name)) {
        auto spd = params.getString(std::string("spectrum " + name));
        return std::make_shared<ConstantTexture const>(
            rgba(Spectrum::loadSPD(pbrtPath + spd)));
    } else if (params.hasString(std::string("texture ") + name)) {
        auto named = params.getString(std::string("texture ") + name);
        auto texture = PBRTState::namedTextures.find(named);
        if (texture != PBRTState::namedTextures.end())
            return texture->second;
        else
            return std::make_shared<ConstantTexture const>(rgba(0, 0, 1, 1));
    } else {
        return std::make_shared<ConstantTexture const>(defaultValue);
    }
}

std::shared_ptr<const Material> makeMaterial(std::string const& type,
                                             Params const& params) {
    std::shared_ptr<const Material> result;

    if (type == "mix") {
        auto m1 = params.getString("string namedmaterial1");
        auto m2 = params.getString("string namedmaterial2");
        auto amount = params.getRGB("rgb amount");

        auto tex = std::make_shared<ConstantTexture const>(amount);
        auto mat1 = PBRTState::namedMaterials.at(m1);
        auto mat2 = PBRTState::namedMaterials.at(m2);

        result = std::make_shared<MixMaterial>(tex, mat1, mat2);
    } else if (type == "glass") {
        auto Kr = makeTexture("Kr", params, rgba(1, 1, 1, 1));
        auto Kt = makeTexture("Kt", params, rgba(1, 1, 1, 1));
        result = std::make_shared<GlassMaterial>(Kr, Kt);
    } else if (type == "metal") {
        printf("Parsing metal\n");
        auto k = makeTexture("k", params, rgba(1, 1, 1, 1));
        auto eta = makeTexture("eta", params, rgba(1, 1, 1, 1));
        result = std::make_shared<MetalMaterial>(eta, k);
    } else {
        result = std::make_shared<DiffuseMaterial>(
            makeTexture("Kd", params, rgba(1, 0, 0, 1)));
    }

    if (params.hasString("texture bumpmap")) {
        auto bumpMap = makeTexture("bumpmap", params, rgba(0, 0, 0, 0));
        result = std::make_shared<BumpMaterial>(bumpMap, result);
    }

    return result;
}

void pbrtMaterial(std::string const& type, Params const& params) {
    PBRTState::materialStack.top() = makeMaterial(type, params);
}

void pbrtMakeNamedMaterial(std::string const& name, Params const& params) {
    auto type = params.getString("string type");
    PBRTState::namedMaterials[name] = makeMaterial(type, params);
}

void pbrtNamedMaterial(std::string const& name) {
    PBRTState::materialStack.top() = PBRTState::namedMaterials.at(name);
}

void pbrtTexture(std::string const& name, std::string const& type,
                 std::string const& clas, Params const& params) {
    if (clas == "imagemap") {
        auto filename = params.getString("string filename");
        std::shared_ptr<const Texture> texture;
        auto cached = PBRTState::cachedTextures.find(pbrtPath + filename);
        if (cached != PBRTState::cachedTextures.end()) {
            texture = cached->second;
        } else {
            texture = std::make_shared<ImageTexture>(pbrtPath + filename);
            PBRTState::cachedTextures[pbrtPath + filename] = texture;
        }

        PBRTState::namedTextures[name] = texture;
    } else if (clas == "scale") {
        std::shared_ptr<Texture const> tex1 =
            makeTexture("tex1", params, rgba(0, 1, 0, 1));
        std::shared_ptr<Texture const> tex2 =
            makeTexture("tex2", params, rgba(0, 1, 0, 1));

        auto texture = std::make_shared<ScaleTexture>(tex1, tex2);
        PBRTState::namedTextures[name] = texture;
    } else {
        std::cout << "Unknown texture class " << clas << std::endl;
    }
}

void pbrtLightSource(std::string const& name, Params const& params) {
    if (name == "infinite") {
        rgba L(1, 1, 1, 1);
        if (params.hasRGB("rgb L")) {
            L = params.getRGB("rgb L");
        }

        std::string mapname;
        if (params.hasString("string mapname")) {
            mapname = pbrtPath + params.getString("string mapname");
        }

        auto light = std::make_shared<InfiniteAreaLight const>(L, mapname);
        PBRTState::objectStack.top().second.lights.push_back(
            LightInstance{light, PBRTState::transforms.top().start});
    }
}

void pbrtWorldBegin() {
    pbrtAttributeBegin();
    pbrtIdentity();
    PBRTState::objectStack.push({"", {}});
}

void pbrtWorldEnd() {
    auto object = PBRTState::objectStack.top();
    PBRTState::objectStack.pop();
    scene =
        std::make_shared<Scene>(object.second.instances, object.second.lights);
    std::cout << "Made World with " << object.second.instances.size()
              << " shapes and " << object.second.lights.size() << " lights"
              << std::endl;
    pbrtAttributeEnd();
}

class Tokens {
   public:
    Tokens(std::vector<std::string> tokens) : tokens(tokens), i(0) {}

    bool done() const { return i >= tokens.size(); }
    std::string peek() const {
        if (!done())
            return tokens.at(i);
        else
            throw;
    }
    std::string next() {
        if (!done())
            return tokens.at(i++);
        else
            throw;
    }

   private:
    std::vector<std::string> tokens;
    size_t i;
};

typedef void (*CommandFunc)(Tokens& tokens);

std::string parseString(Tokens& tokens) {
    auto t = tokens.next();

    // strip off quotes and return
    return t.substr(1, t.length() - 2);
}

float parseReal(Tokens& tokens) { return stof(tokens.next()); }

int64_t parseInt(Tokens& tokens) { return stoll(tokens.next()); }

bool parseBool(Tokens& tokens) {
    auto t = parseString(tokens);
    return t == "true";
}

std::vector<std::string> parseStrings(Tokens& tokens) {
    std::vector<std::string> r;
    if (tokens.peek() == "[") {
        tokens.next();
        while (tokens.peek() != "]") r.push_back(parseString(tokens));
        tokens.next();
    } else {
        r.push_back(parseString(tokens));
    }
    return r;
}

std::vector<float> parseReals(Tokens& tokens) {
    std::vector<float> r;
    if (tokens.peek() == "[") {
        tokens.next();
        while (tokens.peek() != "]") r.push_back(parseReal(tokens));
        tokens.next();
    } else {
        r.push_back(parseReal(tokens));
    }
    return r;
}

std::vector<int64_t> parseInts(Tokens& tokens) {
    std::vector<int64_t> r;
    if (tokens.peek() == "[") {
        tokens.next();
        while (tokens.peek() != "]") r.push_back(parseInt(tokens));
        tokens.next();
    } else {
        r.push_back(parseInt(tokens));
    }
    return r;
}

std::vector<bool> parseBools(Tokens& tokens) {
    std::vector<bool> r;
    if (tokens.peek() == "[") {
        tokens.next();
        while (tokens.peek() != "]") r.push_back(parseBool(tokens));
        tokens.next();
    } else {
        r.push_back(parseBool(tokens));
    }
    return r;
}

Params parseParams(Tokens& tokens) {
    Params r;
    while (tokens.peek()[0] == '"') {
        auto param = parseString(tokens);
        if (param.rfind("string", 0) == 0 || param.rfind("texture", 0) == 0 ||
            param.rfind("spectrum", 0) == 0)
            r.strings[param] = parseStrings(tokens);
        else if (param.rfind("integer", 0) == 0)
            r.ints[param] = parseInts(tokens);
        else if (param.rfind("string", 0) == 0)
            r.bools[param] = parseBools(tokens);
        else
            r.reals[param] = parseReals(tokens);
    }

    return r;
}

void parseCamera(Tokens& tokens) {
    auto name = parseString(tokens);
    auto params = parseParams(tokens);
    pbrtCamera(name, params);
}

void parseSampler(Tokens& tokens) {
    parseString(tokens);
    parseParams(tokens);
}

void parseFilm(Tokens& tokens) {
    parseString(tokens);
    parseParams(tokens);
}

void parsePixelFilter(Tokens& tokens) {
    parseString(tokens);
    parseParams(tokens);
}

void parseIntegrator(Tokens& tokens) {
    parseString(tokens);
    parseParams(tokens);
}

void parseAccelerator(Tokens& tokens) {
    parseString(tokens);
    parseParams(tokens);
}

void parseMakeNamedMedium(Tokens& tokens) {
    parseString(tokens);
    parseParams(tokens);
}

void parseMediumInterface(Tokens& tokens) {
    parseString(tokens);
    parseString(tokens);
}

void parseAttributeBegin(Tokens& tokens) { pbrtAttributeBegin(); }

void parseAttributeEnd(Tokens& tokens) { pbrtAttributeEnd(); }

void parseReverseOrientation(Tokens& tokens) {}

void parseShape(Tokens& tokens) {
    auto name = parseString(tokens);
    auto params = parseParams(tokens);
    pbrtShape(name, params);
}

void parseObjectBegin(Tokens& tokens) {
    auto name = parseString(tokens);
    pbrtObjectBegin(name);
}

void parseObjectEnd(Tokens& tokens) { pbrtObjectEnd(); }

void parseObjectInstance(Tokens& tokens) {
    auto name = parseString(tokens);
    pbrtObjectInstance(name);
}

void parseLightSource(Tokens& tokens) {
    auto name = parseString(tokens);
    auto params = parseParams(tokens);
    pbrtLightSource(name, params);
}

void parseAreaLightSource(Tokens& tokens) {
    parseString(tokens);
    parseParams(tokens);
}

void parseMaterial(Tokens& tokens) {
    auto name = parseString(tokens);
    auto params = parseParams(tokens);
    pbrtMaterial(name, params);
}

void parseMakeNamedMaterial(Tokens& tokens) {
    auto name = parseString(tokens);
    auto params = parseParams(tokens);
    pbrtMakeNamedMaterial(name, params);
}

void parseNamedMaterial(Tokens& tokens) {
    auto name = parseString(tokens);
    pbrtNamedMaterial(name);
}

void parseTexture(Tokens& tokens) {
    auto name = parseString(tokens);
    auto type = parseString(tokens);
    auto clas = parseString(tokens);
    auto params = parseParams(tokens);
    pbrtTexture(name, type, clas, params);
}

void parseWorldBegin(Tokens& tokens) { pbrtWorldBegin(); }

void parseWorldEnd(Tokens& tokens) { pbrtWorldEnd(); }

void parseIdentity(Tokens& tokens) { pbrtIdentity(); }

void parseTranslate(Tokens& tokens) {
    float x = parseReal(tokens);
    float y = parseReal(tokens);
    float z = parseReal(tokens);
    pbrtTranslate(x, y, z);
}

void parseScale(Tokens& tokens) {
    float x = parseReal(tokens);
    float y = parseReal(tokens);
    float z = parseReal(tokens);
    pbrtScale(x, y, z);
}

void parseRotate(Tokens& tokens) {
    float a = parseReal(tokens);
    float x = parseReal(tokens);
    float y = parseReal(tokens);
    float z = parseReal(tokens);
    pbrtRotate(a, x, y, z);
}

void parseLookAt(Tokens& tokens) {
    float ex = parseReal(tokens);
    float ey = parseReal(tokens);
    float ez = parseReal(tokens);
    float lx = parseReal(tokens);
    float ly = parseReal(tokens);
    float lz = parseReal(tokens);
    float ux = parseReal(tokens);
    float uy = parseReal(tokens);
    float uz = parseReal(tokens);
    pbrtLookAt(ex, ey, ez, lx, ly, lz, ux, uy, uz);
}

void parseCoordinateSystem(Tokens& tokens) { parseString(tokens); }

void parseCoordSysTransform(Tokens& tokens) { parseString(tokens); }

void parseTransform(Tokens& tokens) {
    float t[16];
    for (size_t i = 0; i < 16; ++i) t[i] = parseReal(tokens);
    pbrtTransform(t);
}

void parseConcatTransform(Tokens& tokens) {
    float t[16];
    for (size_t i = 0; i < 16; ++i) t[i] = parseReal(tokens);
    pbrtConcatTransform(t);
}

void parseTransformBegin(Tokens& tokens) { pbrtTransformBegin(); }

void parseTransformEnd(Tokens& tokens) { pbrtTransformEnd(); }

void parseActiveTransform(Tokens& tokens) {
    auto s = parseString(tokens);
    if (s == "StartTime")
        pbrtActiveTransformStartTime();
    else if (s == "EndTime")
        pbrtActiveTransformEndTime();
    else if (s == "All")
        pbrtActiveTransformAll();
}

void parseTransformTimes(Tokens& tokens) {
    auto s = parseReal(tokens);
    auto e = parseReal(tokens);
    pbrtTransformTimes(s, e);
}

std::map<std::string, CommandFunc> Commands{
    {"Camera", parseCamera},
    {"Sampler", parseSampler},
    {"Film", parseFilm},
    {"PixelFilter", parsePixelFilter},
    {"Integrator", parseIntegrator},
    {"Accelerator", parseAccelerator},
    {"MakeNamedMedium", parseMakeNamedMedium},
    {"MediumInterface", parseMediumInterface},

    {"AttributeBegin", parseAttributeBegin},
    {"AttributeEnd", parseAttributeEnd},
    {"ReverseOrientation", parseReverseOrientation},
    {"Shape", parseShape},
    {"ObjectBegin", parseObjectBegin},
    {"ObjectEnd", parseObjectEnd},
    {"ObjectInstance", parseObjectInstance},
    {"LightSource", parseLightSource},
    {"AreaLightSource", parseAreaLightSource},
    {"Material", parseMaterial},
    {"MakeNamedMaterial", parseMakeNamedMaterial},
    {"NamedMaterial", parseNamedMaterial},
    {"Texture", parseTexture},

    {"WorldBegin", parseWorldBegin},
    {"WorldEnd", parseWorldEnd},

    {"Identity", parseIdentity},
    {"Translate", parseTranslate},
    {"Scale", parseScale},
    {"Rotate", parseRotate},
    {"LookAt", parseLookAt},
    {"CoordinateSystem", parseCoordinateSystem},
    {"CoordSysTransform", parseCoordSysTransform},
    {"Transform", parseTransform},
    {"ConcatTransform", parseConcatTransform},
    {"TransformBegin", parseTransformBegin},
    {"TransformEnd", parseTransformEnd},
    {"ActiveTransform", parseActiveTransform},
    {"TransformTimes", parseTransformTimes},
};

void parse(Tokens& tokens) {
    while (!tokens.done()) {
        auto token = tokens.next();
        auto c = Commands.find(token);
        if (c != Commands.end())
            c->second(tokens);
        else {
            std::cout << "Unable to find PBRT command: " << token << std::endl;
        }
    }
}

bool iswhitespace(char c) {
    return c == '\n' || c == '\r' || c == '\t' || c == ' ';
}

void tokenize(std::string const& line, std::vector<std::string>& tokens) {
    std::string::size_type i = 0;
    while (i < line.length()) {
        while (i < line.length() && iswhitespace(line[i])) ++i;

        if (i < line.length() && line[i] == '#') i = line.length();

        auto start = i;

        if (i < line.length()) {
            if (line[i] == '[') {
                ++i;
            } else if (line[i] == ']') {
                ++i;
            } else if (line[i] == '"') {
                ++i;
                while (i < line.length() && line[i] != '"') ++i;
                ++i;
            } else if (line[i] == '\'') {
                ++i;
                while (i < line.length() && line[i] != '\'') ++i;
                ++i;
            } else {
                while (i < line.length() && !iswhitespace(line[i]) &&
                       line[i] != ']' && line[i] != '[' && line[i] != '"' &&
                       line[i] != '\'')
                    ++i;
            }

            if (start < i) tokens.push_back(line.substr(start, i - start));
        }
    }
}

void readFile(std::string const& path, std::string const& name,
              std::vector<std::string>& tokens) {
    // std::cout << "Parsing " << path << name << std::endl;

    std::ifstream file(path + name, std::ifstream::in);
    std::string line;
    while (std::getline(file, line)) {
        if (line.rfind("Include", 0) == 0) {
            auto start = line.find("\"", 8);
            auto end = line.find("\"", start + 1);
            if (end != std::string::npos)
                readFile(path, line.substr(start + 1, end - (start + 1)),
                         tokens);
            else
                printf("Error: Include without file name");
        } else {
            tokenize(line, tokens);
        }
    }
}

void parsePBRT(std::string const& file) {
    std::string path = "";
    std::string name = file;
    auto split = file.rfind("/");
    if (split != std::string::npos) {
        path = file.substr(0, split + 1);
        name = file.substr(split + 1);
    }

    std::vector<std::string> tokens;
    readFile(path, name, tokens);
    Tokens tt(tokens);

    pbrtInit(path);
    parse(tt);
    pbrtCleanup();
}
