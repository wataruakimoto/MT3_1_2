#include <Novice.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>

struct Vector3 {
	float x;
	float y;
	float z;
};

struct Matrix4x4 {
	float m[4][4];
};

struct Sphere {
	Vector3 center; // 中心点
	float radius; // 半径
};

// 行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {

	Matrix4x4 resultMultiply = {};

	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			resultMultiply.m[y][x] = m1.m[y][0] * m2.m[0][x] + m1.m[y][1] * m2.m[1][x] + m1.m[y][2] * m2.m[2][x] + m1.m[y][3] * m2.m[3][x];
		}
	}

	return resultMultiply;
}

// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {

	// w = 1 がデカルト座標系であるので(x,y,z,1)のベクトルとしてmatrixとの積をとる
	Vector3 resultTransform = {};

	resultTransform.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
	resultTransform.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
	resultTransform.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];

	// ベクトルに対して基本的な操作を行う行列でwが0になることはありえない
	assert(w != 0.0f);

	// w = 1 がデカルト座標系であるので、w除算することで同次座標をデカルト座標に戻す
	resultTransform.x /= w;
	resultTransform.y /= w;
	resultTransform.z /= w;

	return resultTransform;
}

// 3次元回転行列
Matrix4x4 MakeRotateMatrix(Vector3 radian) {

	Matrix4x4 rotateX = { 0.0f };

	rotateX.m[0][0] = 1.0f;
	rotateX.m[1][1] = cosf(radian.x);
	rotateX.m[1][2] = sinf(radian.x);
	rotateX.m[2][1] = -sinf(radian.x);
	rotateX.m[2][2] = cosf(radian.x);
	rotateX.m[3][3] = 1.0f;

	Matrix4x4 rotateY = { 0.0f };

	rotateY.m[0][0] = cosf(radian.y);
	rotateY.m[0][2] = -sinf(radian.y);
	rotateY.m[1][1] = 1.0f;
	rotateY.m[2][0] = sinf(radian.y);
	rotateY.m[2][2] = cosf(radian.y);
	rotateY.m[3][3] = 1.0f;

	Matrix4x4 rotateZ = { 0.0f };

	rotateZ.m[0][0] = cosf(radian.z);
	rotateZ.m[0][1] = sinf(radian.z);
	rotateZ.m[1][0] = -sinf(radian.z);
	rotateZ.m[1][1] = cosf(radian.z);
	rotateZ.m[2][2] = 1.0f;
	rotateZ.m[3][3] = 1.0f;

	Matrix4x4 resultRotate = { 0.0f };

	resultRotate = Multiply(rotateX, Multiply(rotateY, rotateZ));

	return resultRotate;
}

// 3次元アフィン変換行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {

	Matrix4x4 rotateMatrix = { 0.0f };

	rotateMatrix = MakeRotateMatrix(rotate);

	Matrix4x4 resultAffine = { 0.0f };

	resultAffine.m[0][0] = scale.x * rotateMatrix.m[0][0];
	resultAffine.m[0][1] = scale.x * rotateMatrix.m[0][1];
	resultAffine.m[0][2] = scale.x * rotateMatrix.m[0][2];
	resultAffine.m[1][0] = scale.y * rotateMatrix.m[1][0];
	resultAffine.m[1][1] = scale.y * rotateMatrix.m[1][1];
	resultAffine.m[1][2] = scale.y * rotateMatrix.m[1][2];
	resultAffine.m[2][0] = scale.z * rotateMatrix.m[2][0];
	resultAffine.m[2][1] = scale.z * rotateMatrix.m[2][1];
	resultAffine.m[2][2] = scale.z * rotateMatrix.m[2][2];

	resultAffine.m[3][0] = translate.x;
	resultAffine.m[3][1] = translate.y;
	resultAffine.m[3][2] = translate.z;
	resultAffine.m[3][3] = 1.0f;

	return resultAffine;
}

// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {

	Matrix4x4 resultPerspectiveFov = {};

	resultPerspectiveFov.m[0][0] = (1 / aspectRatio) * (1 / tanf(fovY / 2));
	resultPerspectiveFov.m[1][1] = 1 / tanf(fovY / 2);
	resultPerspectiveFov.m[2][2] = farClip / (farClip - nearClip);
	resultPerspectiveFov.m[2][3] = 1.0f;
	resultPerspectiveFov.m[3][2] = -nearClip * farClip / (farClip - nearClip);

	return resultPerspectiveFov;
}

// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {

	Matrix4x4 resultPerspectiveFov = {};

	resultPerspectiveFov.m[0][0] = (1 / aspectRatio) * (1 / tanf(fovY / 2));
	resultPerspectiveFov.m[1][1] = 1 / tanf(fovY / 2);
	resultPerspectiveFov.m[2][2] = farClip / (farClip - nearClip);
	resultPerspectiveFov.m[2][3] = 1.0f;
	resultPerspectiveFov.m[3][2] = -nearClip * farClip / (farClip - nearClip);

	return resultPerspectiveFov;
}

// ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {

	Matrix4x4 resultViewport = {};

	resultViewport.m[0][0] = width / 2;
	resultViewport.m[1][1] = -height / 2;
	resultViewport.m[2][2] = maxDepth - minDepth;
	resultViewport.m[3][0] = left + (width / 2);
	resultViewport.m[3][1] = top + (height / 2);
	resultViewport.m[3][2] = minDepth;
	resultViewport.m[3][3] = 1.0f;

	return resultViewport;
}

// 球の描画
void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {

	const uint32_t kSubdivision = 6; // 分割数
	const float kLonEvery = 180 / kSubdivision; // 経度分割1つ分の角度 30
	const float kLatEvery = 180 / kSubdivision; // 緯度分割1つ分の角度 30

	// 緯度の方向に分割 -Π/2 ~ Π/2
	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {

		float lat = -M_PI / 2.0f + kLatEvery * latIndex; // 現在の緯度

		// 経度の方向に分割 0 ~ 2Π
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {

			float lon = lonIndex * kLonEvery; // 現在の経度

			// World座標系でのa、b、cを求める
			Vector3 a, b, c;

			a = {
				(sphere.center.x + sphere.radius) * cosf(lat) * cosf(lon),
				(sphere.center.y + sphere.radius) * sinf(lat),
				(sphere.center.z + sphere.radius) * cosf(lat) * sinf(lon)
			};

			b = {
				(sphere.center.x + sphere.radius) * cosf(lat + (M_PI / kSubdivision)) * cosf(lon),
				(sphere.center.y + sphere.radius) * sinf(lat + (M_PI / kSubdivision)),
				(sphere.center.z + sphere.radius) * cosf(lat + (M_PI / kSubdivision)) * sinf(lon)
			};

			c = {
				(sphere.center.x + sphere.radius) * cosf(lat) * cosf(lon + ((M_PI * 2) / kSubdivision)),
				(sphere.center.y + sphere.radius) * sinf(lat),
				(sphere.center.z + sphere.radius) * cosf(lat) * sinf(lon + ((M_PI * 2) / kSubdivision))
			};

			// a、b、cをScreen座標系まで変換
			Matrix4x4 worldAMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, a);
			Matrix4x4 worldBMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, b);
			Matrix4x4 worldCMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, c);

			Matrix4x4 worldViewProjectionAMatrix = Multiply(worldAMatrix, viewProjectionMatrix);
			Matrix4x4 worldViewProjectionBMatrix = Multiply(worldAMatrix, viewProjectionMatrix);
			Matrix4x4 worldViewProjectionCMatrix = Multiply(worldAMatrix, viewProjectionMatrix);

			Vector3 ndcAMatrix = Transform(a, worldViewProjectionAMatrix);
			Vector3 ndcBMatrix = Transform(b, worldViewProjectionBMatrix);
			Vector3 ndcCMatrix = Transform(c, worldViewProjectionCMatrix);

			// ab、bcで線を引く
			Novice::DrawLine()
		}
	}
}

// グリッドの描画
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {

	const float kGridHalfWidth = 2.0f;                                      // Gridの半分の幅
	const uint32_t kSubdivision = 10;                                       // 分割数
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision); // 1つ分の長さ

	// 奥から手前への線を順々にに引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {

	}
}

const char kWindowTitle[] = "LE2B_01_アキモト_ワタル";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
