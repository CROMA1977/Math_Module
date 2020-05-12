#include <math.h>
#include <vector>
#include <algorithm>

#include "CVectorReference.h"
#include "CVector.h"

// 建構式
CVectorReference2FLOAT::CVectorReference2FLOAT(float & Value0, float & Value1) :
m_x(Value0), m_y(Value1)
{
}

// 建構式 Array[2] = {m_x, m_y}
CVectorReference2FLOAT::CVectorReference2FLOAT(float * Array) :
m_x(Array[0]), m_y(Array[1])
{
}

CVectorReference2FLOAT::CVectorReference2FLOAT(CVector2FLOAT & Vector) :
m_x(Vector.m_x), m_y(Vector.m_y)
{
}

// 拷貝建構式
CVectorReference2FLOAT::CVectorReference2FLOAT(const CVectorReference2FLOAT & Vector) :
m_x(Vector.m_x), m_y(Vector.m_y)
{
}

// 解構式
CVectorReference2FLOAT::~CVectorReference2FLOAT()
{
}

// 設定向量 
CVectorReference2FLOAT & CVectorReference2FLOAT::Set(const float x, const float y)
{
	m_x = x; m_y = y;
	return *this;
}

// 設定向量 Array[2] = {m_x, m_y}
CVectorReference2FLOAT & CVectorReference2FLOAT::Set(const float * Array)
{
	m_x = Array[0]; m_y = Array[1];
	return *this;
}

// 設定向量 V = V0 + (V1-V0) * t 
CVectorReference2FLOAT & CVectorReference2FLOAT::Set(const CVectorReference2FLOAT & V0, const CVectorReference2FLOAT & V1, const float t)
{
	// V0 + (V1-V0) * t 
	// V0 + V1 * t - V0 * t
	// V0 - V0 * t + V1 * t
	// V0 * (1-t) + V1 * t
	const float one_minus_t = 1 - t;
	m_x = V0.m_x * one_minus_t + V1.m_x * t;
	m_y = V0.m_y * one_minus_t + V1.m_y * t;
	return *this;
}

// 計算兩向量的 Dot (Ax * Bx + Ay * By)
float CVectorReference2FLOAT::Dot(const CVectorReference2FLOAT & Vector) const
{
	return m_x * Vector.m_x + m_y * Vector.m_y;
}

// 計算兩向量的方向 逆時針 > 0 順時針 < 0
float CVectorReference2FLOAT::Cross(const CVectorReference2FLOAT & Vector) const
{
	return m_x * Vector.m_y - m_y * Vector.m_x;
}

// 使兩向量交換 
void CVectorReference2FLOAT::Swap(CVectorReference2FLOAT & Vector)
{
	CVectorReference2FLOAT Temp(Vector);
	Vector = *this;
	*this = Temp;
}

// 取得兩向量的弳度
float CVectorReference2FLOAT::GetRadian(const CVectorReference2FLOAT & Vector) const
{
	float Length1 = (*this).GetLength();
	float Length2 = Vector.GetLength();
	float Length = Length1 * Length2;
	if (Length == 0) {
		fprintf(stderr, "\n error:Vec2f GetAngle Length == 0!");
		return 0;
	}
	// 限制 在 -1 ~ 1 之間
	float Cos = max(-1.0f, min(1.0f, (*this).Dot(Vector) / Length));
	return acosf(Cos);
}

// 取得兩向量的角度
float CVectorReference2FLOAT::GetAngle(const CVectorReference2FLOAT & Vector) const
{
	return GetRadian(Vector) * 180.0f / 3.1415926f;
}

// 取得向量的長度
float CVectorReference2FLOAT::GetLength(void) const
{
	return sqrtf(m_x * m_x + m_y * m_y);
}

// 取得向量長度的平方
float CVectorReference2FLOAT::GetLengthSquare(void) const
{
	return m_x * m_x + m_y * m_y;
}

// 取得單位化的向量
CVector2FLOAT CVectorReference2FLOAT::GetNormalize(void) const
{
	float Length = GetLength();
	if (Length == 0) {
		fprintf(stderr, "\n error:Vec2f GetNormalize Length == 0!");
		return CVector2FLOAT(*this);
	}
	return CVector2FLOAT(m_x / Length, m_y / Length);
}

// 取得反轉的向量 V[2] = {-m_x, -m_y}
CVector2FLOAT CVectorReference2FLOAT::GetReverse(void) const
{
	return CVector2FLOAT(-m_x, -m_y);
}

// 使向量單位化
CVectorReference2FLOAT & CVectorReference2FLOAT::Normalize(void)
{
	float Length = GetLength();
	if (Length == 0) {
		fprintf(stderr, "\n error:Vec2f Normalize Length == 0!");
		return *this;
	}
	m_x /= Length;
	m_y /= Length;
	return *this;
}
// 使向量反轉 
CVectorReference2FLOAT& CVectorReference2FLOAT::Reverse(void)
{
	m_x = -m_x;	m_y = -m_y;
	return *this;
}

// 檢查向量是否為零
bool CVectorReference2FLOAT::IsZero(void) const
{
	float LengthPow2 = m_x * m_x + m_y * m_y;
	if (LengthPow2 > 0.000001)
		return false;
	return true;
}
// 覆載對向量的 + 運算子 V[2] = [m_x + Vector.m_x, m_y + Vector.m_y]
CVector2FLOAT CVectorReference2FLOAT::operator+(const CVectorReference2FLOAT& Vector) const
{
	return CVector2FLOAT( m_x + Vector.m_x, m_y + Vector.m_y );
}

// 覆載對向量的 - 運算子 V[2] = [m_x - Vector.m_x, m_y - Vector.m_y]
CVector2FLOAT CVectorReference2FLOAT::operator-(const CVectorReference2FLOAT& Vector) const
{
	return CVector2FLOAT(m_x - Vector.m_x, m_y - Vector.m_y);
}

// 覆載對實數的 * 運算子 V[2] = [m_x * Value, m_y * Value]
CVector2FLOAT CVectorReference2FLOAT::operator*(const float Value) const
{
	return CVector2FLOAT(m_x * Value, m_y * Value);
}

// 覆載對實數的 / 運算子 V[2] = [m_x / Value, m_y / Value] 
CVector2FLOAT CVectorReference2FLOAT::operator/(const float Value) const
{
	if (Value == 0) {
		fprintf(stderr, "\n error:Vec2f Div Value == 0!");
		return CVector2FLOAT(*this);
	}
	return CVector2FLOAT(m_x / Value, m_y / Value);
}

// 覆載 = 運算子 V = Vector
CVectorReference2FLOAT & CVectorReference2FLOAT::operator=(const CVectorReference2FLOAT& Vector)
{
	m_x = Vector.m_x; m_y = Vector.m_y;
	return *this;
}
// 覆載 = 運算子 V = Vector
CVectorReference2FLOAT & CVectorReference2FLOAT::operator=(const CVector2FLOAT & Vector)
{
	m_x = Vector.m_x; m_y = Vector.m_y;
	return *this;
}

// 覆載對向量的 += 運算子 [m_x, m_y] = [m_x + Vector.m_x, m_y + Vector.m_y]
CVectorReference2FLOAT& CVectorReference2FLOAT::operator+=(const CVectorReference2FLOAT& Vector)
{
	m_x += Vector.m_x; m_y += Vector.m_y;
	return *this;
}

// 覆載對向量的 -= 運算子 [m_x, m_y] = [m_x - Vector.m_x, m_y - Vector.m_y] 
CVectorReference2FLOAT& CVectorReference2FLOAT::operator-=(const CVectorReference2FLOAT& Vector)
{
	m_x -= Vector.m_x; m_y -= Vector.m_y;
	return *this;
}

// 覆載對實數的 *= 運算子 [m_x, m_y] = [m_x * Value, m_y * Value]
CVectorReference2FLOAT& CVectorReference2FLOAT::operator*=(const float Value)
{
	m_x *= Value; m_y *= Value;
	return *this;
}

// 覆載對實數的 /= 運算子 [m_x, m_y] = [m_x / Value, m_y / Value]
CVectorReference2FLOAT & CVectorReference2FLOAT::operator/=(const float Value)
{
	if (Value == 0) {
		fprintf(stderr, "\n error:Vec2f DivAssign Value == 0!");
		return *this;
	}
	m_x /= Value;
	m_y /= Value;
	return *this;
}

// 覆載對向量的 == 運算子 dx * dx + dy * dy <  0.000001
bool CVectorReference2FLOAT::operator==(const CVectorReference2FLOAT& Vector) const
{
	float dx = m_x - Vector.m_x;
	float dy = m_y - Vector.m_y;
	float LengthPow2 = dx * dx + dy * dy;

	if (LengthPow2 >= 0.000001)
		return false;
	return true;
}

bool CVectorReference2FLOAT::operator!=(const CVectorReference2FLOAT& Vector) const
{
	float dx = m_x - Vector.m_x;
	float dy = m_y - Vector.m_y;
	float LengthPow2 = dx * dx + dy * dy;

	if (LengthPow2 >= 0.000001)
		return true;
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

// 建構式
CVectorReference3FLOAT::CVectorReference3FLOAT(float & Value0, float & Value1, float & Value2) :
m_x(Value0), m_y(Value1), m_z(Value2)
{

}

// 建構式 Array[3] = {m_x, m_y, m_z};
CVectorReference3FLOAT::CVectorReference3FLOAT(float * Array) :
m_x(Array[0]), m_y(Array[1]), m_z(Array[2])
{

}

CVectorReference3FLOAT::CVectorReference3FLOAT(CVector3FLOAT& Vector) :
m_x(Vector.m_x), m_y(Vector.m_y), m_z(Vector.m_z)
{
}

// 拷貝建構式	
CVectorReference3FLOAT::CVectorReference3FLOAT(const CVectorReference3FLOAT & Vector) :
m_x(Vector.m_x), m_y(Vector.m_y), m_z(Vector.m_z)
{

}

// 解構式
CVectorReference3FLOAT::~CVectorReference3FLOAT() {

}

// 設定向量
CVectorReference3FLOAT & CVectorReference3FLOAT::Set(const float x, const float y, const float z) {
	m_x = x; m_y = y; m_z = z;
	return *this;
}

// 設定向量 Array[3] = {m_x, m_y, m_z}
CVectorReference3FLOAT & CVectorReference3FLOAT::Set(const float * Array) {
	m_x = Array[0]; m_y = Array[1]; m_z = Array[2];
	return *this;
}

// 設定向量 V = V0 + (V1-V0) * t 
CVectorReference3FLOAT& CVectorReference3FLOAT::Set(const CVectorReference3FLOAT& V0, const CVectorReference3FLOAT& V1, const float t) {
	// V0 + (V1-V0) * t 
	// V0 + V1 * t - V0 * t
	// V0 - V0 * t + V1 * t
	// V0 * (1-t) + V1 * t
	const float one_minus_t = 1 - t;
	m_x = V0.m_x * one_minus_t + V1.m_x * t;
	m_y = V0.m_y * one_minus_t + V1.m_y * t;
	m_z = V0.m_z * one_minus_t + V1.m_z * t;
	return *this;
}

// 計算兩向量的 Dot (Ax * Bx + Ay * By + Az * Bz)
float CVectorReference3FLOAT::Dot(const CVectorReference3FLOAT& Vector) const {
	return m_x * Vector.m_x + m_y * Vector.m_y + m_z * Vector.m_z;
}

// 計算兩向量的 Cross 
CVector3FLOAT CVectorReference3FLOAT::Cross(const CVectorReference3FLOAT & Vector) const {
	float X, Y, Z;

	// 計算結果
	X = m_y * Vector.m_z - m_z * Vector.m_y;
	Y = m_z * Vector.m_x - m_x * Vector.m_z;
	Z = m_x * Vector.m_y - m_y * Vector.m_x;

	return CVector3FLOAT(X, Y, Z);
}

// 使兩向量交換 
void CVectorReference3FLOAT::Swap(CVectorReference3FLOAT& Vector) {
	CVectorReference3FLOAT Temp(Vector);
	Vector = *this;
	*this = Temp;
}

// 取得兩向量的弳度
float CVectorReference3FLOAT::GetRadian(const CVectorReference3FLOAT & Vector) const {
	float LengthA = (*this).GetLength();
	float LengthB = Vector.GetLength();
	float Length = LengthA * LengthB;
	// 任何一個向量的長度為零時回傳角度為 0
	if (Length == 0) {
		fprintf(stderr, "\n error:Vec3f GetAngle Length == 0!");
		return 0;
	}
	// ( A Dot B ) / (La * Lb) 限制 在 -1 ~ 1 之間
	float Cos = max(-1.0f, min(1.0f, (*this).Dot(Vector) / Length));
	// 轉換成弳度
	return acosf(Cos);
}

// 取得兩向量的角度
float CVectorReference3FLOAT::GetAngle(const CVectorReference3FLOAT & Vector) const {
	// 轉換成角度
	return GetRadian(Vector) * 180 / 3.1415926f;
}

// 取得向量的長度
float CVectorReference3FLOAT::GetLength(void) const {
	float LengthPow2 = m_x * m_x + m_y * m_y + m_z * m_z;
	return sqrtf(LengthPow2);
}

// 取得向量長度的平方
float CVectorReference3FLOAT::GetLengthSquare(void) const {
	return m_x * m_x + m_y * m_y + m_z * m_z;
}

// 取得單位化的向量
CVector3FLOAT CVectorReference3FLOAT::GetNormalize(void) const {
	float Length = GetLength();
	if (Length == 0) {
		fprintf(stderr, "\n error:Vec3f GetNormalize Length == 0!");
		return CVector3FLOAT(0.0f, 0.0f, 0.0f);
	}
	return CVector3FLOAT( m_x / Length, m_y / Length, m_z / Length );
}

// 取得反轉的向量 V[3] = {-m_x, -m_y, -m_z}
CVector3FLOAT CVectorReference3FLOAT::GetReverse(void) const {
	return CVector3FLOAT(-m_x, -m_y, -m_z);
}

// 使向量單位化
CVectorReference3FLOAT & CVectorReference3FLOAT::Normalize(void) {
	float Length = GetLength();
	if (Length == 0) {
		fprintf(stderr, "\n error:Vec3f Normalize Length == 0!");
		return *this;
	}
	m_x /= Length; m_y /= Length; m_z /= Length;
	return *this;
}

// 使向量反轉
CVectorReference3FLOAT& CVectorReference3FLOAT::Reverse(void) {
	m_x = -m_x; m_y = -m_y; m_z = -m_z;
	return *this;
}

// 檢查向量是否為零
bool CVectorReference3FLOAT::IsZero(void) const {
	double LengthPow2 = m_x * m_x + m_y * m_y + m_z * m_z;
	if (LengthPow2 > 0.000001)
		return false;
	return true;
}

// 覆載對向量的 + 運算子 V[3] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z]
CVector3FLOAT CVectorReference3FLOAT::operator+(const CVectorReference3FLOAT& Vector) const {
	return CVector3FLOAT( m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z );
}

// 覆載對向量的 - 運算子 V[3] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z]
CVector3FLOAT CVectorReference3FLOAT::operator-(const CVectorReference3FLOAT& Vector) const {
	return CVector3FLOAT( m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z );
}

// 覆載對實數的 * 運算子 V[3] = [m_x * Value, m_y * Value, m_z * Value]
CVector3FLOAT CVectorReference3FLOAT::operator*(const float Value) const {
	return CVector3FLOAT( m_x * Value, m_y * Value, m_z * Value );
}

// 覆載對實數的 / 運算子 V[3] = [m_x / Value, m_y / Value, m_z / Value] 
CVector3FLOAT CVectorReference3FLOAT::operator/(const float Value) const {
	if (Value == 0) {
		fprintf(stderr, "\n error:Vec3f Div Value == 0!");
		return CVector3FLOAT(*this);
	}
	return CVector3FLOAT( m_x / Value, m_y / Value, m_z / Value);
}

// 覆載 = 運算子 V = Vector
CVectorReference3FLOAT& CVectorReference3FLOAT::operator=(const CVectorReference3FLOAT & Vector) {
	m_x = Vector.m_x; m_y = Vector.m_y; m_z = Vector.m_z;
	return *this;
}

// 覆載 = 運算子 V = Vector
CVectorReference3FLOAT & CVectorReference3FLOAT::operator=(const CVector3FLOAT & Vector)
{
	m_x = Vector.m_x; m_y = Vector.m_y; m_z = Vector.m_z;
	return *this;
}

// 覆載對向量的 += 運算子 [m_x, m_y, m_z] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z]
CVectorReference3FLOAT& CVectorReference3FLOAT::operator+=(const CVectorReference3FLOAT & Vector) {
	m_x += Vector.m_x; m_y += Vector.m_y; m_z += Vector.m_z;
	return *this;
}

// 覆載對向量的 -= 運算子 [m_x, m_y, m_z] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z] 
CVectorReference3FLOAT& CVectorReference3FLOAT::operator-=(const CVectorReference3FLOAT & Vector) {
	m_x -= Vector.m_x; m_y -= Vector.m_y; m_z -= Vector.m_z;
	return *this;
}

// 覆載對實數的 *= 運算子 [m_x, m_y, m_z] = [m_x * Value, m_y * Value, m_z * Value]
CVectorReference3FLOAT& CVectorReference3FLOAT::operator*=(const float Value) {
	m_x *= Value; m_y *= Value; m_z *= Value;
	return *this;
}

// 覆載對實數的 /= 運算子 [m_x, m_y, m_z] = [m_x / Value, m_y / Value, m_z / Value] 
CVectorReference3FLOAT& CVectorReference3FLOAT::operator/=(const float Value) {
	if (Value == 0) {
		fprintf(stderr, "\n error:Vec3f DivAssign Value == 0!");
		return *this;
	}
	m_x /= Value; m_y /= Value; m_z /= Value;
	return *this;
}

// 覆載對向量的 == 運算子 dx * dx + dy * dy + dz * dz < 0.000001
bool CVectorReference3FLOAT::operator==(const CVectorReference3FLOAT& Vector) const {
	float dx = m_x - Vector.m_x;
	float dy = m_y - Vector.m_y;
	float dz = m_z - Vector.m_z;
	float LengthPow2 = dx * dx + dy * dy + dz * dz;

	if (LengthPow2 >= 0.000001)
		return false;
	return true;
}

// 覆載對向量的 != 運算子 dx * dx + dy * dy + dz * dz >= 0.000001
bool CVectorReference3FLOAT::operator!=(const CVectorReference3FLOAT& Vector) const {
	float dx = m_x - Vector.m_x;
	float dy = m_y - Vector.m_y;
	float dz = m_z - Vector.m_z;
	float LengthPow2 = dx * dx + dy * dy + dz * dz;

	if (LengthPow2 >= 0.000001)
		return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////

// 建構式
CVectorReference4FLOAT::CVectorReference4FLOAT(float & Value0, float & Value1, float & Value2, float & Value3) :
m_x(Value0), m_y(Value1), m_z(Value2), m_w(Value3)
{
}

// 建構式
CVectorReference4FLOAT::CVectorReference4FLOAT(float * Array) :
m_x(Array[0]), m_y(Array[1]), m_z(Array[2]), m_w(Array[3])
{
}

CVectorReference4FLOAT::CVectorReference4FLOAT(CVector4FLOAT & Vector) :
m_x(Vector.m_x), m_y(Vector.m_y), m_z(Vector.m_z), m_w(Vector.m_w)
{
}

// 拷貝建構式 
CVectorReference4FLOAT::CVectorReference4FLOAT(const CVectorReference4FLOAT & Vector) :
m_x(Vector.m_x), m_y(Vector.m_y), m_z(Vector.m_z), m_w(Vector.m_w)
{
}

// 解構式
CVectorReference4FLOAT::~CVectorReference4FLOAT() {
}

// 設定向量
CVectorReference4FLOAT& CVectorReference4FLOAT::Set(const float x, const float y, const float z, const float w) {
	m_x = x; m_y = y; m_z = z; m_w = w;
	return *this;
}

// 設定向量 Array[4] = {m_x, m_y, m_z, m_w}
CVectorReference4FLOAT& CVectorReference4FLOAT::Set(const float * Array) {
	m_x = Array[0]; m_y = Array[1]; m_z = Array[2]; m_w = Array[3];
	return *this;
}

// 使兩向量交換 
void CVectorReference4FLOAT::Swap(CVectorReference4FLOAT & Vector) {
	// 複製 Vector 的值
	CVectorReference4FLOAT Temp(Vector);
	Vector = *this;
	*this = Temp;
}

// 取得反轉的向量 V[4] = {-m_x, -m_y, -m_z, -m_w}
CVector4FLOAT CVectorReference4FLOAT::GetReverse(void) const {
	return CVector4FLOAT(-m_x, -m_y, -m_z, -m_w);
}

// 使向量反轉
CVectorReference4FLOAT & CVectorReference4FLOAT::Reverse(void) {
	m_x = -m_x; m_y = -m_y; m_z = -m_z; m_w = -m_w;
	return *this;
}

// 檢查向量是否為零
bool CVectorReference4FLOAT::IsZero(void) const {
	float LengthPow2 = m_x * m_x + m_y * m_y + m_z * m_z + m_w * m_w;
	if (LengthPow2 > 0.000001)
		return false;
	return true;
}

// 覆載對向量的 + 運算子 V[4] = [m_x + Vector.x(), m_y + Vector.y(), m_z + Vector.z(), m_w + Vector.w()]
CVector4FLOAT CVectorReference4FLOAT::operator+(const CVectorReference4FLOAT& Vector) const {
	return CVector4FLOAT(m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z, m_w + Vector.m_w);
}

// 覆載對向量的 - 運算子 V[4] = [m_x - Vector.x(), m_y - Vector.y(), m_z - Vector.z(), m_w - Vector.w()]
CVector4FLOAT CVectorReference4FLOAT::operator-(const CVectorReference4FLOAT& Vector) const {
	return CVector4FLOAT(m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z, m_w - Vector.m_w);
}

// 覆載對實數的 * 運算子 V[4] = [m_x * Value, m_y * Value, m_z * Value, m_w * Value]
CVector4FLOAT CVectorReference4FLOAT::operator*(const float Value) const {
	return CVector4FLOAT(m_x * Value, m_y * Value, m_z * Value, m_w * Value);
}

// 覆載對實數的 / 運算子 V[4] = [m_x / Value, m_y / Value, m_z / Value, m_w / Value] 
CVector4FLOAT CVectorReference4FLOAT::operator/(const float Value) const {
	if (Value == 0) {
		fprintf(stderr, "\n error:Vec4f Div Value == 0!");
		return CVector4FLOAT(*this);
	}
	return CVector4FLOAT(m_x / Value, m_y / Value, m_z / Value, m_w / Value);
}

// 覆載 = 運算子 V = Vector
CVectorReference4FLOAT& CVectorReference4FLOAT::operator=(const CVectorReference4FLOAT & Vector) {
	m_x = Vector.m_x; m_y = Vector.m_y; m_z = Vector.m_z; m_w = Vector.m_w;
	return *this;
}

// 覆載 = 運算子 V = Vector
CVectorReference4FLOAT & CVectorReference4FLOAT::operator=(const CVector4FLOAT & Vector)
{
	m_x = Vector.m_x; m_y = Vector.m_y; m_z = Vector.m_z; m_w = Vector.m_w;
	return *this;
}

// 覆載對向量的 += 運算子 [m_x, m_y, m_z, m_w] = [m_x + Vector.x(), m_y + Vector.y(), m_z + Vector.z(), m_w + Vector.w()]
CVectorReference4FLOAT& CVectorReference4FLOAT::operator+=(const CVectorReference4FLOAT& Vector) {
	m_x += Vector.m_x; m_y += Vector.m_y; m_z += Vector.m_z; m_w += Vector.m_w;
	return *this;
}

// 覆載對向量的 -= 運算子 [m_x, m_y, m_z, m_w] = [m_x - Vector.x(), m_y - Vector.y(), m_z - Vector.z(), m_w - Vector.w()]
CVectorReference4FLOAT& CVectorReference4FLOAT::operator-=(const CVectorReference4FLOAT& Vector) {
	m_x -= Vector.m_x; m_y -= Vector.m_y; m_z -= Vector.m_z; m_w -= Vector.m_w;
	return *this;
}

// 覆載對實數的 *= 運算子 [m_x, m_y, m_z, m_w] = [m_x * Value, m_y * Value, m_z * Value, m_w * Value]
CVectorReference4FLOAT& CVectorReference4FLOAT::operator*=(const float Value) {
	m_x *= Value; m_y *= Value; m_z *= Value; m_w *= Value;
	return *this;
}

// 覆載對實數的 /= 運算子 [m_x, m_y, m_z, m_w] = [m_x / Value, m_y / Value, m_z / Value, m_w / Value] 
CVectorReference4FLOAT& CVectorReference4FLOAT::operator/=(const float Value) {
	if (Value == 0) {
		fprintf(stderr, "\n error:Vec4f DivAssign Value == 0!");
		return *this;
	}
	m_x /= Value; m_y /= Value; m_z /= Value; m_w /= Value;
	return *this;
}

// 覆載對向量的 == 運算子 dx * dx + dy * dy + dz * dz + dw * dw < 0.000001
bool CVectorReference4FLOAT::operator==(const CVectorReference4FLOAT& Vector) const {
	float dx = m_x - Vector.m_x;
	float dy = m_y - Vector.m_y;
	float dz = m_z - Vector.m_z;
	float dw = m_w - Vector.m_w;
	float LengthPow2 = dx * dx + dy * dy + dz * dz + dw * dw;
	if (LengthPow2 >= 0.000001)
		return false;
	return true;
}

// 覆載對向量的 != 運算子 dx * dx + dy * dy + dz * dz + dw * dw >= 0.000001
bool CVectorReference4FLOAT::operator!=(const CVectorReference4FLOAT& Vector) const {
	float dx = m_x - Vector.m_x;
	float dy = m_y - Vector.m_y;
	float dz = m_z - Vector.m_z;
	float dw = m_w - Vector.m_w;
	float LengthPow2 = dx * dx + dy * dy + dz * dz + dw * dw;
	if (LengthPow2 >= 0.000001)
		return true;
	return false;
}


//////////////////////////////////////////////////////////////////////////