export class ResponseCodedError extends Error {
  constructor(public error_code: number, message?: string) {
    super(message)
    Object.setPrototypeOf(this, ResponseCodedError.prototype);
  }
}

/**
 * The method is not supported
 */
export class UnsupportedMethodError extends ResponseCodedError {
  constructor(message = 'Unsupported method') {
    super(405, message)
    Object.setPrototypeOf(this, UnsupportedMethodError.prototype)
  }
}

/**
 * The resource is not found
 */
export class NotFoundError extends ResponseCodedError {
  constructor(message = 'Not Found') {
    super(404, message)
    Object.setPrototypeOf(this, NotFoundError.prototype)
  }
}

/**
 * The resource is not found
 */
export class UnauthorizedError extends ResponseCodedError {
  constructor(message = 'Unauthorized') {
    super(401, message)
    Object.setPrototypeOf(this, UnauthorizedError.prototype)
  }
}

/**
 * This timeout error is used to ensure we don't wait too long for dependencies
 *  fortunately, even if it occurs the job will requeue still making progress.
 */
export class TimeoutError extends ResponseCodedError {
  constructor(message = 'Timeout reached') {
    super(504, message)
    Object.setPrototypeOf(this, TimeoutError.prototype)
  }
}

/**
 * This error occurs when the input node is not populated yet
 */
export class UnboundError extends ResponseCodedError {
  constructor(message = 'Refusing to submit unbound variable') {
    super(422, message)
    Object.setPrototypeOf(this, UnboundError.prototype)
  }
}
