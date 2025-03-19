//! By convention, root.zig is the root source file when making a library. If
//! you are making an executable, the convention is to delete this file and
//! start with main.zig instead.
const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const ArrayList = std.ArrayList;
const ff = std.crypto.ff;
const print = std.debug.print;

pub fn UnivariatePolynomial(comptime M: type, comptime T: type, x: comptime_int) type {
    return struct {
        const Self = @This();

        /// The coefficient of `x^i` is stored at location `i` in `self.coeffs`.
        coefficients: ArrayList(M.Fe),
        m: M,

        /// Initializes a new univariate polynomial from a `slice` of coefficients of primitive type `T`.
        /// The `slice` must have been allocated with `allocator`.
        ///
        /// Deinitialize with `deinit` or use `toOwnedSlice`.
        pub fn fromOwnedSlicePrimitive(allocator: Allocator, slice: []T) Self {
            var coefficients = ArrayList(M.Fe).initCapacity(allocator, slice.len) catch unreachable;
            const m = M.fromPrimitive(T, x) catch unreachable;
            for (slice) |c| {
                const fe = M.Fe.fromPrimitive(T, m, c) catch unreachable;
                coefficients.appendAssumeCapacity(fe);
            }
            return .{
                .coefficients = coefficients,
                .m = m,
            };
        }

        /// Initializes a new univariate polynomial from a `slice` of coefficients of type `M.Fe`.
        /// The `slice` must have been allocated with `allocator`.
        ///
        /// Deinitialize with `deinit` or use `toOwnedSlice`.
        pub fn fromOwnedSlice(allocator: Allocator, slice: []M.Fe) Self {
            const coefficients = ArrayList(M.Fe).fromOwnedSlice(allocator, slice);
            return .{
                .coefficients = coefficients,
                .m = M.fromPrimitive(T, x) catch unreachable,
            };
        }

        /// Returns the degree of the univariate polynomial.
        pub fn degree(self: Self) usize {
            return if (self.coefficients.items.len == 0) 0 else self.coefficients.items.len - 1;
        }

        /// Horner's method to evaluate polynomials
        /// Takes n additions and n multiplications
        fn evaluate(self: Self, point: M.Fe) M.Fe {
            var result = self.coefficients.items[self.coefficients.items.len - 1];
            for (1..self.coefficients.items.len) |i| {
                result = self.m.add(self.m.mul(result, point), self.coefficients.items[self.coefficients.items.len - 1 - i]);
            }

            return result;
        }

        fn deinit(self: Self) void {
            self.coefficients.deinit();
        }
    };
}

test "basic polynomial functionality fromPrimitive" {
    const M = ff.Modulus(256);
    const m = try M.fromPrimitive(u256, 31);
    const expected = try M.Fe.fromPrimitive(u256, m, 17);

    var coefficients: [3]u256 = [_]u256{ 1, 2, 3 };
    const poly = UnivariatePolynomial(M, u256, 31).fromOwnedSlicePrimitive(std.testing.allocator, coefficients[0..]);
    defer poly.deinit();
    const point = try M.Fe.fromPrimitive(u256, m, 2);
    const result = poly.evaluate(point);

    try std.testing.expectEqual(expected, result);
}

test "basic polynomial functionality" {
    const M = ff.Modulus(256);
    const m = try M.fromPrimitive(u256, 31);
    const x = try M.Fe.fromPrimitive(u256, m, 1);
    const y = try M.Fe.fromPrimitive(u256, m, 2);
    const z = try M.Fe.fromPrimitive(u256, m, 3);
    const expected = try M.Fe.fromPrimitive(u256, m, 17);

    var coefficients = try std.testing.allocator.alloc(M.Fe, 3);
    coefficients[0] = x;
    coefficients[1] = y;
    coefficients[2] = z;
    const poly = UnivariatePolynomial(M, u256, 31).fromOwnedSlice(std.testing.allocator, coefficients[0..]);
    defer poly.deinit();
    const point = try M.Fe.fromPrimitive(u256, m, 2);
    const result = poly.evaluate(point);

    try std.testing.expectEqual(expected, result);
}
